/** 
FEM Vector Fields
 */
#ifdef WIN32
#define NOMINMAX
#include <windows.h>
#endif

#if defined (__APPLE__) || defined (OSX)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include "GA/gl_util.h"

#include "primitivedraw.h"
#include "gahelper.h"
#include "Laplacian.h"

#include <memory>

#include <vector>
#include <queue>
#include <map>
#include <fstream>
#include <functional>
#include <complex>
#include "numerics.h"
#include "HalfEdge/Mesh.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

// #include <ppl.h>

const char *WINDOW_TITLE = "FEM BASIC 2D";

// GLUT state information
int g_viewportWidth = 800;
int g_viewportHeight = 600;

void display();
void reshape(GLint width, GLint height);
void MouseButton(int button, int state, int x, int y);
void MouseMotion(int x, int y);
void KeyboardUpFunc(unsigned char key, int x, int y);
void SpecialFunc(int key, int x, int y);
void SpecialUpFunc(int key, int x, int y);
void Idle();
void DestroyWindow();
Eigen::Vector3d valueToColor( double d );

//using namespace boost;
using namespace c3ga;
using namespace std;
using namespace numerics;

class Camera
{
public:
	float		pos[3];
	float		fw[3];
	float		up[3];
	float		translateVel;
	float		rotateVel;

	Camera()
	{
		float		_pos[] = { 0, 0, 2};
		float		_fw[] = { 0, 0, -1 };
		float		_up[] = { 0, 1, 0 };

		translateVel = 0.005;
		rotateVel = 0.005;
		memcpy(pos, _pos, sizeof(float)*3);
		memcpy(fw, _fw, sizeof(float)*3);
		memcpy(up, _up, sizeof(float)*3);
	}

	void glLookAt()
	{
		gluLookAt( pos[0], pos[1], pos[2], fw[0],  fw[1],  fw[2], up[0],  up[1],  up[2] );
	}
};

class VertexBuffer
{
public:
	std::vector<Eigen::Vector3d> positions; //mesh vertex positions
	std::vector<Eigen::Vector3d> normals; //for rendering (lighting)
	std::vector<Eigen::Vector3d> colors; //for rendering (visual representation of values)
	int size;

	VertexBuffer() : size(0)
	{
	}

	void resize(int size)
	{
		this->size = size;
		positions.resize(size);
		normals.resize(size);
		colors.resize(size);
	}
	int get_size() { return size; }

};

class IndexBuffer {
public:
	std::vector<int> faces;
	int size;

	IndexBuffer() : size(0)
	{
	}

	void resize(int size)
	{
		this->size = size;
		faces.resize(size);
	}
	int get_size() { return size; }

};

Camera g_camera;
Mesh mesh;
vectorE3GA g_prevMousePos;
bool g_rotateModel = false;
bool g_rotateModelOutOfPlane = false;
rotor g_modelRotor = _rotor(1.0);
float g_dragDistance = -1.0f;
int g_dragObject;
bool g_showWires = true;


VertexBuffer vertexBuffer;
IndexBuffer indexBuffer;

/**
 * Compute normal per vertex and tangent plane
 * Compute a basis tangent vectors per vertex
 * Compute a normal vector per face
 * Compute a 2x2 matrix aligning face normal with vertex normal
 * Use that matrix to compute the connection-laplacian
 */

typedef Eigen::Matrix<double, 3, 2> Matrix32d;
typedef Eigen::Matrix<double, 2, 3> Matrix23d;
typedef Eigen::Matrix<double, 4, 6> Matrix46d;
typedef Eigen::Matrix<double, 6, 6> Matrix66d;
std::vector<Matrix32d> vertex_tangent_space;
std::vector<double> face_areas;
std::vector<Eigen::Vector3d> face_normals;
std::vector<Matrix32d> face_tangent_space;
std::shared_ptr<SparseMatrix<double>> A;
Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
std::set<int> allconstraints;
Eigen::VectorXd right_hand_side;
Eigen::VectorXd solutionU;

/// Project u on the orthgonal of n
/// \param u vector to project
/// \param n vector to build orthogonal space from
/// \return projected vector
static Eigen::Vector3d project(const Eigen::Vector3d & u, const Eigen::Vector3d & n)
{
	return u - (u.dot(n) / n.squaredNorm()) * n;
}

void computeVertexTangentSpace(Mesh *mesh, vector<Matrix32d>& vertex_tangent_space) {
	Matrix32d tangentSpace;
	for(Vertex &vi : mesh->getVertices()) {
		Eigen::Vector3d& ni = vi.n;
		Eigen::Vector3d& pi = vi.p;
		Eigen::Vector3d& pj = vi.edge->pair->vertex->p;
		Eigen::Vector3d tangent = project(pj - pi, ni).normalized();
		Eigen::Vector3d bitangent = ni.cross(tangent);
		tangentSpace.col(0) = tangent;
		tangentSpace.col(1) = bitangent;
		vertex_tangent_space[vi.ID] = tangentSpace;
	}
}

void computeFaceNormals(Mesh *mesh, vector<Eigen::Vector3d>& face_normals, vector<double>& face_areas) {
	Vertex* v[3];
	for(Face &face : mesh->getFaces()) {
		v[0] = face.edge->vertex;
		v[1] = face.edge->next->vertex;
		v[2] = face.edge->next->next->vertex;
		Eigen::Vector3d uu = v[1]->p - v[0]->p;
		Eigen::Vector3d vv = v[2]->p - v[0]->p;
		Eigen::Vector3d n = uu.cross(vv);
		face_areas[face.ID] = 0.5 * n.norm();
		face_normals[face.ID] = n.normalized();
	}
}

void computeFaceTangentSpace(Mesh *mesh, const vector<Eigen::Vector3d>& face_normals, vector<Matrix32d>& face_tangent_space) {
	Matrix32d tangentSpace;
	Vertex* v[2];
	for(Face &face : mesh->getFaces()) {
		v[0] = face.edge->vertex;
		v[1] = face.edge->next->vertex;
		Eigen::Vector3d uu = (v[1]->p - v[0]->p).normalized();
		tangentSpace.col(0) = uu;
		tangentSpace.col(1) = face_normals[face.ID].cross(uu);
		face_tangent_space[face.ID] = tangentSpace;
	}
}

///Return [n] as the 3x3 operator such that [n]q = n x q
///@param n a vector
Eigen::Matrix3d bracket(const Eigen::Vector3d &n)
{
	Eigen::Matrix3d brack;
	brack << 0.0,  -n(2), n(1),
             n(2), 0.0 , -n(0),
            -n(1), n(0), 0.0 ;
	return brack;
}


/**
 * Computes the gradient of the shape function "N^k_i" at the nodes "i" on a reference 2D triagle "k"
 * then "push forward" the gradients to the pysical triangle using the "differential" of the mapping.
 * All of that follows from the chain-rule as described in:
 * 
 * Francisco-Javier Sayas, "A gentle introduction to the Finite Element Method", 2015
 */
Eigen::Matrix3d Gradient(const Face& f) {

	Matrix32d B;
	Matrix32d Binv;
	Eigen::Matrix3d gradN;
	Vertex *v[3];
	v[0] = f.edge->vertex;
	v[1] = f.edge->next->vertex;
	v[2] = f.edge->next->next->vertex;

	B(0,0) = v[1]->p.x() - v[0]->p.x(); B(0,1) = v[2]->p.x() - v[0]->p.x();
	B(1,0) = v[1]->p.y() - v[0]->p.y(); B(1,1) = v[2]->p.y() - v[0]->p.y();
	B(2,0) = v[1]->p.z() - v[0]->p.z(); B(2,1) = v[2]->p.z() - v[0]->p.z();
    
	Binv = ((B.transpose() * B).inverse() * B.transpose()).transpose();

	//grad N^k_1(X) = B^-T (-1, -1)
	//grad N^k_2(X) = B^-T (1, 0)
	//grad N^k_3(X) = B^-T (0, 1)

	gradN.col(0) = Eigen::Vector3d(-Binv(0,0) - Binv(0,1), -Binv(1,0) - Binv(1,1), -Binv(2,0) - Binv(2,1));
	gradN.col(1) = Eigen::Vector3d( Binv(0,0), Binv(1,0), Binv(2,0));
	gradN.col(2) = Eigen::Vector3d( Binv(0,1), Binv(1,1), Binv(2,1));

	return gradN; 
}


/// https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
///@return 3x3 Rotation matrix to align n_v to n_f
Eigen::Matrix3d alignVectors(const Eigen::Vector3d& nv, const Eigen::Vector3d& nf)
{
	double c = nv.dot(nf);
	//Special case for opposite nv and nf vectors.
	if (std::abs(c + 1.0) < 0.00001)
		return -Eigen::Matrix3d::Identity();

	auto vv = nv.cross(nf);
	Eigen::Matrix3d skew = bracket(vv);
	return Eigen::Matrix3d::Identity() + skew + 1.0 / (1.0 + c) * skew * skew;
}


///@return Levi-Civita connection from vertex v tangent space to face f
///tangent space (2x2 rotation matrix)
Eigen::Matrix2d Rvf(const Vertex* v, const Face* f)
{
	return face_tangent_space[f->ID].transpose() * alignVectors(v->n, face_normals[f->ID]) * vertex_tangent_space[v->ID];
}

///@return Levi-Civita connection from vertex v tangent space to face f
///tangent space (2x2 rotation matrix)
Eigen::Matrix2d Rvv(const Vertex* vi, const Vertex* vj)
{
	return vertex_tangent_space[vj->ID].transpose() * alignVectors(vi->n, vj->n) * vertex_tangent_space[vi->ID];
}

/**
 * Extend computation of per-element stiffness matrix of triangle elements embedded in 3D space.
 * Original method only works for triangle elements on 2D space, see:
 * 
 * Francisco-Javier Sayas, "A gentle introduction to the Finite Element Method", 2015
 */
Matrix66d AssembleStiffnessElementEmbedded(Vertex* v[3]) {

	Eigen::MatrixXd B(3, 2);
	Eigen::MatrixXd Binv(3, 2);
	Eigen::Vector3d gradN[3];
	Matrix66d elementMatrix;

	B(0,0) = v[1]->p.x() - v[0]->p.x(); B(0,1) = v[2]->p.x() - v[0]->p.x();
	B(1,0) = v[1]->p.y() - v[0]->p.y(); B(1,1) = v[2]->p.y() - v[0]->p.y();
	B(2,0) = v[1]->p.z() - v[0]->p.z(); B(2,1) = v[2]->p.z() - v[0]->p.z();
    
	Binv = ((B.transpose() * B).inverse() * B.transpose()).transpose();

	double faceArea = 0.5 * ((v[1]->p - v[0]->p).cross(v[2]->p - v[0]->p)).norm();

	//grad N^k_1(X) = B^-T (-1, -1)
	//grad N^k_2(X) = B^-T (1, 0)
	//grad N^k_3(X) = B^-T (0, 1)

	gradN[0] = Eigen::Vector3d(-Binv(0,0) - Binv(0,1), -Binv(1,0) - Binv(1,1), -Binv(2,0) - Binv(2,1));
	gradN[1] = Eigen::Vector3d( Binv(0,0), Binv(1,0), Binv(2,0));
	gradN[2] = Eigen::Vector3d( Binv(0,1), Binv(1,1), Binv(2,1));
	elementMatrix.setZero();
	for( int i = 0 ; i < 3 ; ++i ) { // for each test function
		for (int j = 0 ; j < 3 ; ++j ) { // for each shape function
			if (i < j) continue; // since stifness matrix is symmetric
			//w_ij = area K <grad N^k_i(X), grad N^k_j(X)>
			if(i == j) {
				elementMatrix.block<2,2>(i * 2, j * 2) = faceArea * gradN[i].dot(gradN[j]) * Eigen::Matrix2d::Identity();
			}
			else {
				elementMatrix.block<2,2>(i * 2, j * 2) = faceArea * gradN[i].dot(gradN[j]) * Rvv(v[j], v[i]);
				elementMatrix.block<2,2>(j * 2, i * 2) = faceArea * gradN[j].dot(gradN[i]) * Rvv(v[i], v[j]);
			}
		}
	}
	return elementMatrix;
}

/**
 * Extend computation of per-element stiffness matrix of triangle elements embedded in 3D space.
 * Original method only works for triangle elements on 2D space, see:
 * 
 * Francisco-Javier Sayas, "A gentle introduction to the Finite Element Method", 2015
 */
Matrix66d AssembleStiffnessElementEmbedded2(Vertex* v[3], const Face& face) {

	Eigen::MatrixXd B(3, 2);
	Eigen::MatrixXd Binv(3, 2);
	Eigen::Vector3d gradN[3];
	Matrix66d elementMatrix;

	B(0,0) = v[1]->p.x() - v[0]->p.x(); B(0,1) = v[2]->p.x() - v[0]->p.x();
	B(1,0) = v[1]->p.y() - v[0]->p.y(); B(1,1) = v[2]->p.y() - v[0]->p.y();
	B(2,0) = v[1]->p.z() - v[0]->p.z(); B(2,1) = v[2]->p.z() - v[0]->p.z();
    
	Binv = ((B.transpose() * B).inverse() * B.transpose()).transpose();

	double faceArea = 0.5 * ((v[1]->p - v[0]->p).cross(v[2]->p - v[0]->p)).norm();

	//grad N^k_1(X) = B^-T (-1, -1)
	//grad N^k_2(X) = B^-T (1, 0)
	//grad N^k_3(X) = B^-T (0, 1)

	gradN[0] = Eigen::Vector3d(-Binv(0,0) - Binv(0,1), -Binv(1,0) - Binv(1,1), -Binv(2,0) - Binv(2,1));
	gradN[1] = Eigen::Vector3d( Binv(0,0), Binv(1,0), Binv(2,0));
	gradN[2] = Eigen::Vector3d( Binv(0,1), Binv(1,1), Binv(2,1));
	elementMatrix.setZero();

	Matrix32d Tf = face_tangent_space[face.ID];

	for( int i = 0 ; i < 3 ; ++i ) { // for each test function
		for (int j = 0 ; j < 3 ; ++j ) { // for each shape function
			if (i < j) continue; // since stifness matrix is symmetric
			//w_ij = area K <grad N^k_i(X), grad N^k_j(X)>
			Eigen::Vector2d grad_i = Tf.transpose() * gradN[i];
			Eigen::Vector2d grad_j = Tf.transpose() * gradN[j];
			Eigen::Matrix2d R_i = Rvf(v[i], &face);
			Eigen::Matrix2d R_j = Rvf(v[j], &face);
			
			if(i == j) {
				elementMatrix.block<2,2>(i * 2, j * 2) = faceArea * grad_i.dot(grad_j) * Eigen::Matrix2d::Identity();
				
			}
			else {
				elementMatrix.block<2,2>(i * 2, j * 2) = faceArea * grad_i.dot(grad_j) * R_i.transpose() * R_j;
				elementMatrix.block<2,2>(j * 2, i * 2) = faceArea * grad_j.dot(grad_i) * R_j.transpose() * R_i;
			}
		}
	}
	return elementMatrix;
}

/**
 * Assemble the stiffness matrix. It does not take into account boundary conditions.
 * Boundary conditions will be applied when linear system is pre-factored (LU decomposition)
 * Original method can be found in:
 * 
 * Francisco-Javier Sayas, "A gentle introduction to the Finite Element Method", 2015
 */
std::shared_ptr<SparseMatrix<double>> AssembleMatrix(Mesh *mesh, double delta_t) {
	std::shared_ptr<SparseMatrix<double>> A(new SparseMatrix<double>(mesh->numVertices() * 2, mesh->numVertices() * 2));
	Matrix66d stiffnessMatrix;
	Eigen::Matrix2d wij;
	Vertex* v[3];
	for (Face& face : mesh->getFaces()) {
		v[0] = face.edge->vertex;
		v[1] = face.edge->next->vertex;
		v[2] = face.edge->next->next->vertex;
		stiffnessMatrix = AssembleStiffnessElementEmbedded(v);
		for( int i = 0 ; i < 3 ; ++i ) {
			for (int j = 0 ; j < 3 ; ++j ) {
				//wij = massMatrix(i, j) + delta_t * stiffnessMatrix(i, j);
				wij = delta_t * stiffnessMatrix.block<2,2>(i * 2, j * 2);
				for(int row = 0 ; row < 2 ; ++row) {
					for(int col = 0 ; col < 2 ; ++col) {
						//if(wij(row, col) != 0.0) 
						{
							(*A)(v[i]->ID * 2 + row, v[j]->ID * 2 + col) += wij(row, col);
						}
					}
				}
			}
		}
	}
	return A;
}

std::shared_ptr<SparseMatrix<double>> AssembleMatrix2(Mesh *mesh, double delta_t) {
	std::shared_ptr<SparseMatrix<double>> A(new SparseMatrix<double>(mesh->numVertices() * 2, mesh->numVertices() * 2));
	Matrix66d stiffnessMatrix;
	Eigen::Matrix2d wij;
	Vertex* v[3];
	for (Face& face : mesh->getFaces()) {
		v[0] = face.edge->vertex;
		v[1] = face.edge->next->vertex;
		v[2] = face.edge->next->next->vertex;
		//stiffnessMatrix = connectionLaplacian(face);
		stiffnessMatrix = AssembleStiffnessElementEmbedded2(v, face);
		for( int i = 0 ; i < 3 ; ++i ) {
			for (int j = 0 ; j < 3 ; ++j ) {
				//wij = massMatrix(i, j) + delta_t * stiffnessMatrix(i, j);
				wij = delta_t * stiffnessMatrix.block<2,2>(i * 2, j * 2);
				for(int row = 0 ; row < 2 ; ++row) {
					for(int col = 0 ; col < 2 ; ++col) {
						//if(wij(row, col) != 0.0) 
						{
							(*A)(v[i]->ID * 2 + row, v[j]->ID * 2 + col) += wij(row, col);
						}
					}
				}
			}
		}
	}
	return A;
}



double computeEdgeLengths(
	Mesh* mesh
) {
	int count = 0;
	double sumEdgeLengths = 0.0;
	for(auto &eij : mesh->getEdges()) {
		Eigen::Vector3d& pi = eij->vertex->p;
		Eigen::Vector3d& pj = eij->pair->vertex->p;
		sumEdgeLengths += (pj - pi).norm();
		count++;
	}
	return sumEdgeLengths / (double)count;
}

std::shared_ptr<SparseMatrix<double>> AssembleDiagonalMassMatrix(Mesh *mesh) {
	std::shared_ptr<SparseMatrix<double>> A(new SparseMatrix<double>(mesh->numVertices() * 2, mesh->numVertices() * 2));
	double wij;
	Vertex* v[3];
	for (Face& face : mesh->getFaces()) {
		v[0] = face.edge->vertex;
		v[1] = face.edge->next->vertex;
		v[2] = face.edge->next->next->vertex;
		double faceArea = face_areas[face.ID];
		wij = faceArea / 3.0;
		(*A)(v[0]->ID * 2    , v[0]->ID * 2    ) += wij;
		(*A)(v[0]->ID * 2 + 1, v[0]->ID * 2 + 1) += wij;
		(*A)(v[1]->ID * 2    , v[1]->ID * 2    ) += wij;
		(*A)(v[1]->ID * 2 + 1, v[1]->ID * 2 + 1) += wij;
		(*A)(v[2]->ID * 2    , v[2]->ID * 2    ) += wij;
		(*A)(v[2]->ID * 2 + 1, v[2]->ID * 2 + 1) += wij;
	}
	return A;
}

void IplusMinvTimesA(std::shared_ptr<SparseMatrix<double>> M, std::shared_ptr<SparseMatrix<double>> A)
{
	auto numRows = A->numRows();
	for (int i = 0; i < numRows; ++i)
	{
		SparseMatrix<double>::RowIterator aIter = A->iterator(i);
		double oneOverVertexOneRingArea = 1.0 / (*M)(i, i);
		for (; !aIter.end(); ++aIter)
		{
			auto j = aIter.columnIndex();
			(*A)(i, j) *= oneOverVertexOneRingArea;
			if (i == j) {
				(*A)(i, j) += 1.0; // this completes the (I + M^-1 L)
			}
		}
	}
}

bool is_constrained(std::set<int>& constraints, int vertex)
{
	return constraints.find(vertex) != constraints.end();
}

void PreFactor(std::shared_ptr<SparseMatrix<double>> A, std::set<int>& constraints, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
{

	Eigen::SparseMatrix<double> Lc = Eigen::SparseMatrix<double>(A->numRows(), A->numColumns());

	auto numRows = A->numRows();
	for (int i = 0; i < numRows; ++i)
	{
		if (!is_constrained(constraints, i))
		{
			SparseMatrix<double>::RowIterator aIter = A->iterator(i);
			for (; !aIter.end(); ++aIter)
			{
				auto j = aIter.columnIndex();
				Lc.insert(i, j) = (*A)(i, j);
			}
		}
		else
		{
			Lc.insert(i, i) = 1.0;
		}
	}

	Lc.makeCompressed();
	solver.compute(Lc);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: " << "Prefactor failed." << std::endl;
		exit(1);
	}
}

int main(int argc, char* argv[])
{
	/**
	 * Load the FEM mesh
	 */
	//mesh.readFEM("lake_nodes.txt", "lake_elements.txt");
	mesh.readOBJ("cactus1.obj");
	mesh.CenterAndNormalize();
	mesh.computeNormals();

	// GLUT Window Initialization:
	glutInit (&argc, argv);
	glutInitWindowSize(g_viewportWidth, g_viewportHeight);
	glutInitDisplayMode( GLUT_RGB | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow(WINDOW_TITLE);

	// Register callbacks:
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(MouseButton);
	glutMotionFunc(MouseMotion);
	glutKeyboardUpFunc(KeyboardUpFunc);
	glutSpecialFunc(SpecialFunc);
	glutSpecialUpFunc(SpecialUpFunc);
	glutIdleFunc(Idle);
	atexit(DestroyWindow);

	InitializeDrawing();

	vertexBuffer.resize(mesh.numVertices());
	indexBuffer.resize(mesh.numFaces() * 3);

	/**
	 * Initialize the vertex-buffer for OpenGL rendering purposes
	 */
	for( Vertex& vertex : mesh.getVertices())
	{
		vertexBuffer.positions[vertex.ID] = vertex.p;
		vertexBuffer.normals[vertex.ID] = vertex.n;
		vertexBuffer.colors[vertex.ID] = valueToColor(0);
	}

	/**
	 * Initialize the index-buffer for OpenGL rendering purposes
	 */
	for (Face& face : mesh.getFaces()) {
		int i = face.ID;
		int	v1 = face.edge->vertex->ID;
		int	v2 = face.edge->next->vertex->ID;
		int	v3 = face.edge->next->next->vertex->ID;
		indexBuffer.faces[i * 3 + 0] = v1;
		indexBuffer.faces[i * 3 + 1] = v2;
		indexBuffer.faces[i * 3 + 2] = v3;
	}

	vertex_tangent_space.resize(mesh.numVertices());
	face_areas.resize(mesh.numFaces());
	face_normals.resize(mesh.numFaces());
	face_tangent_space.resize(mesh.numFaces());

	computeVertexTangentSpace(&mesh, vertex_tangent_space);
	computeFaceNormals(&mesh, face_normals, face_areas);
	computeFaceTangentSpace(&mesh, face_normals, face_tangent_space);

	double avgEdgeLength = computeEdgeLengths(&mesh);

	right_hand_side = Eigen::VectorXd(mesh.numVertices() * 2);
	right_hand_side.setZero(); // solve laplace's equation where RHS is zero

	for( Vertex& vertex : mesh.getVertices())
	{
		if (vertex.p.norm() < 2.5e-2) {
			std::complex<double> U_i = std::complex<double>(cos(2*M_PI/3), sin(2*M_PI/3));
			right_hand_side(vertex.ID * 2    ) = U_i.real();
			right_hand_side(vertex.ID * 2 + 1) = U_i.imag();
			allconstraints.insert(vertex.ID * 2    );
			allconstraints.insert(vertex.ID * 2 + 1);
			break;
		}
		if (vertex.p.z() > 0.83) {
			std::complex<double> U_i = std::complex<double>(cos(2*M_PI/3), sin(2*M_PI/3));
			right_hand_side(vertex.ID * 2    ) = U_i.real();
			right_hand_side(vertex.ID * 2 + 1) = U_i.imag();
			allconstraints.insert(vertex.ID * 2    );
			allconstraints.insert(vertex.ID * 2 + 1);
			break;
		}
	}

	A = AssembleMatrix2(&mesh, avgEdgeLength*avgEdgeLength);
	std::shared_ptr<SparseMatrix<double>> M;
	M = AssembleDiagonalMassMatrix(&mesh);
	IplusMinvTimesA(M, A);

	PreFactor(A, allconstraints, solver);

	solutionU = solver.solve(right_hand_side);

	glutMainLoop();

	return 0;
}

void display()
{
	/*
	 *	matrices
	 */
	glViewport( 0, 0, g_viewportWidth, g_viewportHeight );
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	pickLoadMatrix();
	GLpick::g_frustumFar = 1000.0;
	GLpick::g_frustumNear = .1;
	gluPerspective( 60.0, (double)g_viewportWidth/(double)g_viewportHeight, GLpick::g_frustumNear, GLpick::g_frustumFar );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glShadeModel(GL_SMOOTH);	//gouraud shading
	glClearDepth(1.0f);
	glClearColor( .75f, .75f, .75f, .0f );
	glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

	/*
	 *	estados
	 */
	glEnable(GL_CULL_FACE);		//face culling
	glCullFace( GL_BACK );
	glFrontFace( GL_CCW );
	glEnable(GL_DEPTH_TEST);	//z-buffer
	glDepthFunc(GL_LEQUAL);

	/*
	 *	iluminacion
	 */
	float		ambient[] = { .3f, .3f, .3f, 1.f };
	float		diffuse[] = { .3f, .3f, .3f, 1.f };
	float		position[] = { .0f, 0.f, 15.f, 1.f };
	float		specular[] = { 1.f, 1.f, 1.f };

	glLightfv( GL_LIGHT0, GL_AMBIENT, ambient );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse );
	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0125);
	glEnable(  GL_LIGHT0   );
	glEnable(  GL_LIGHTING );
	//glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specular );
	glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 50.f );

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glLoadIdentity();

	g_camera.glLookAt();

	glLightfv( GL_LIGHT0, /*GL_SPOT_DIRECTION*/GL_POSITION, position );

	glPushMatrix();

	rotorGLMult(g_modelRotor);

	if (GLpick::g_pickActive) glLoadName((GLuint)-1);

	double alpha = 1.0;

	//glEnable (GL_BLEND);
	//glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//alpha = 0.5;

	//Mesh-Faces Rendering
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL /*GL_LINE GL_FILL GL_POINT*/);
	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable( GL_COLOR_MATERIAL );
	if (GLpick::g_pickActive) glLoadName((GLuint)10);

	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, &vertexBuffer.positions[0]);
	glNormalPointer(GL_DOUBLE, 0, &vertexBuffer.normals[0]);
	glColorPointer(3, GL_DOUBLE, 0, &vertexBuffer.colors[0]);

	// draw the model
	glDrawElements(GL_TRIANGLES, indexBuffer.get_size(), GL_UNSIGNED_INT, &indexBuffer.faces[0]);
	// deactivate vertex arrays after drawing
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	if (g_showWires)
	{
		if (!GLpick::g_pickActive)
		{
			//Mesh-Edges Rendering (superimposed to faces)
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE /*GL_LINE GL_FILL GL_POINT*/);
			glColor4d(.5, .5, .5, alpha);
			glDisable(GL_LIGHTING);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_DOUBLE, 0, &vertexBuffer.positions[0]);
			// draw the model
			glDrawElements(GL_TRIANGLES, indexBuffer.get_size(), GL_UNSIGNED_INT, &indexBuffer.faces[0]);
			// deactivate vertex arrays after drawing
			glDisableClientState(GL_VERTEX_ARRAY);
			glEnable(GL_LIGHTING);
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL /*GL_LINE GL_FILL GL_POINT*/);

		}
	}

	glDisable( GL_COLOR_MATERIAL );
	glDisable(GL_POLYGON_OFFSET_FILL);


	//glDisable (GL_BLEND);
	float		green[] = { .0f, .5f, .0f, 1.f };
	float		red[] = { .5f, .0f, .0f, 1.f };
	float		blue[] = { .0f, .0f, .5f, 1.f };

	for(int i = 0; i < mesh.numVertices(); ++i){
		Vertex& vi = mesh.vertexAt(i);
		std::complex<double> R_i = std::complex<double>(solutionU[vi.ID * 2], solutionU[vi.ID * 2 + 1]);
		R_i = R_i / abs(R_i);
		Eigen::Vector3d X_i = project(vi.edge->pair->vertex->p - vi.p, vi.n).normalized();
		Eigen::Quaterniond Q(Eigen::AngleAxisd(log(R_i).imag(), vi.n));
		Eigen::Vector3d U_i = 0.06 * Q._transformVector(X_i);
		DrawArrow(c3gaPoint(vi.p.x(), vi.p.y(), vi.p.z()), _vectorE3GA(U_i.x(), U_i.y(), U_i.z()));
	}

	glPopMatrix();

	glutSwapBuffers();
}

Eigen::Vector3d valueToColor( double d )
{
	static Eigen::Vector3d	c0 = Eigen::Vector3d( 1, 1, 1);
	static Eigen::Vector3d	c1 = Eigen::Vector3d( 1, 1, 0);
	static Eigen::Vector3d	c2 = Eigen::Vector3d( 0, 1, 0);
	static Eigen::Vector3d	c3 = Eigen::Vector3d( 0, 1, 1);
	static Eigen::Vector3d	c4 = Eigen::Vector3d( 0, 0, 1);

	if( d < 0.25 )
	{
		double alpha = (d - 0.0) / (0.25-0.0);
		return (1.0 - alpha) * c0 + alpha * c1;
	}
	else if( d < 0.5 )
	{
		double alpha = (d - 0.25) / (0.5-0.25);
		return (1.0 - alpha) * c1 + alpha * c2;
	}
	else if( d < 0.75 )
	{
		double alpha = (d - 0.5) / (0.75-0.5);
		return (1.0 - alpha) * c2 + alpha * c3;
	}
	else
	{
		double alpha = (d - 0.75) / (1.0-0.75);
		return (1.0 - alpha) * c3 + alpha * c4;
	}
}


void reshape(GLint width, GLint height)
{
	g_viewportWidth = width;
	g_viewportHeight = height;

	// redraw viewport
	glutPostRedisplay();
}

vectorE3GA mousePosToVector(int x, int y) {
	x -= g_viewportWidth / 2;
	y -= g_viewportHeight / 2;
	return _vectorE3GA((float)-x * e1 - (float)y * e2);
}

void MouseButton(int button, int state, int x, int y)
{
	g_rotateModel = false;

	if (button == GLUT_LEFT_BUTTON)
	{
		g_prevMousePos = mousePosToVector(x, y);

		GLpick::g_pickWinSize = 1;
		g_dragObject = pick(x, g_viewportHeight - y, display, &g_dragDistance);

		if(g_dragObject == -1 || g_dragObject == 10 )
		{
			vectorE3GA mousePos = mousePosToVector(x, y);
			g_rotateModel = true;

			if ((_Float(norm_e(mousePos)) / _Float(norm_e(g_viewportWidth * e1 + g_viewportHeight * e2))) < 0.2)
				g_rotateModelOutOfPlane = true;
			else g_rotateModelOutOfPlane = false;
		}
	}

	if (button == GLUT_RIGHT_BUTTON)
	{
		g_prevMousePos = mousePosToVector(x, y);

		GLpick::g_pickWinSize = 1;
		g_dragObject = pick(x, g_viewportHeight - y, display, &g_dragDistance);
	}
}

void MouseMotion(int x, int y)
{
	if (g_rotateModel )
	{
		// get mouse position, motion
		vectorE3GA mousePos = mousePosToVector(x, y);
		vectorE3GA motion = mousePos - g_prevMousePos;

		if (g_rotateModel)
		{
			// update rotor
			if (g_rotateModelOutOfPlane)
				g_modelRotor = exp(g_camera.rotateVel * (motion ^ e3) ) * g_modelRotor;
			else 
				g_modelRotor = exp(0.00001f * (motion ^ mousePos) ) * g_modelRotor;
		}

		// remember mouse pos for next motion:
		g_prevMousePos = mousePos;

		// redraw viewport
		glutPostRedisplay();
	}
}

void SpecialFunc(int key, int x, int y)
{
	switch(key) {
		case GLUT_KEY_F1 :
			{
				int mod = glutGetModifiers();
				if(mod == GLUT_ACTIVE_CTRL || mod == GLUT_ACTIVE_SHIFT )
				{
				}
			}
			break;
		case GLUT_KEY_UP:
			{
			}
			break;
		case GLUT_KEY_DOWN:
			{
			}
			break;
	}
}

void SpecialUpFunc(int key, int x, int y)
{
}

void KeyboardUpFunc(unsigned char key, int x, int y)
{
	if(key == 'w' || key == 'W')
	{
		g_showWires = !g_showWires;
		glutPostRedisplay();
	}
}

void Idle()
{
	// redraw viewport
}

void DestroyWindow()
{
	ReleaseDrawing();
}

