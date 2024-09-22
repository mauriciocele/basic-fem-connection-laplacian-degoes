/** 
FEM Vector Fields


FEM


Solve Poisson equation:

- ConnectionLaplacian U(x,y) = F(x,y)   in the region D
                                        U(x,y) = G(x,y)  on the region boundary #D

 with Dirichlet boundary conditions (boundary value problem)

The ConnectionLaplacian is defined in terms of the “covariant derivative” Nabla as:

ConnectionLaplacian( U(x,y) ) = Nabla^* Nabla ( U(x,y) )

Where Nabla^* is the Adjoint of Nabla (covariant derivative), U(x,y) is a tangent vector field on the 3D surface 
(2D manifold embedded in R^3). F(x,y) is a vector field resulting of applying the ConnectionLaplacian and G(x,y) 
the vector boundary conditions.

Nabla(U, w) is an operator that has two arguments: U(x,y) and w(x,y). U(x,y) is a vector field which is differentiated 
and w(x,y) is a vector field telling the direction in which directional derivatives are taken. And returns a vector field as a 
result which is a measure of relative change of U w.r.t to w. 

Usually one can omit the vector field w and just talk about Nabla U which is a 2-tensor (i.e., a matrix). For example in flat 
space R^3 the vector field represented as U = [u_x, u_y, u_z]^T the operator Nabla U would be:


Nabla U =  d/dx    u_x  u_y  u_z   =   d/dx u_x    d/dy u_x      d/dz u_x
           d/dy                        d/dx u_y    d/dy u_y      d/dz u_y
           d/dz                        d/dx u_z    d/dy u_z      d/dz u_z

Nabla U = Grad Operator  @  U = Jacobian (U)

Where “@“ is the symbol I use for “tensor product“

In a curved manifold this definition would need to be adjusted to cope with curvature with the aid of a “connection”.

The “scalar product” of rank-2 tensors X and Y is denoted as X : Y and is defined as 

X : Y = Trace(X^T Y)

The Adjoint of Nabla is defined by “integration by parts” (see wikipedia https://en.wikipedia.org/wiki/Vector_calculus_identities) 
over surface D:

Integral_D Nabla U : X dx = Integral_D <U, Nabla^* X> dx + Integral_#D <U, X(N)> dA    (1)

where N is unit normal at boundary #D.

The solution U(x,y) should also satisfy the integral form (notice Integral_D is a double integral):

- Integral_D( ConnectionLaplacian U(x,y) dA ) = Integral_D ( F(x,y) dA )

where ConnectionLaplacian U(x,y) =  Nabla^* Nabla ( U(x,y) )


WEAK FORM


Given a test function V(x,y) with compact support which is a vector field, we multiply the PDE by V(x,y) and integrate:

-  Integral_D <Nabla^* Nabla ( U(x,y) ), V(x,y)> dA = Integral_D <F(x,y), V(x,y)> dA                                   (2)

Replacing X with Nabla U in (1) we get:

Integral_D Nabla(V) : Nabla(U) dx = Integral_D <V, Nabla^* Nabla(U)> dx + Integral_#D <V, Nabla(U, N)> dA    (3)

where N is unit normal at boundary #D.

Replacing (3) in (2)

- Integral_D Nabla(V) : Nabla(U) dx = Integral_D <F(x,y), V(x,y)> dA  - Integral_#D <V, Nabla(U, N)> dA

We now define that all V(x, y) at the boundary #D has zero value, so

V(x,y) = 0 for all (x,y) in #D

Which gives:

- Integral_D Nabla(V) : Nabla(U) dx = Integral_D <F(x,y), V(x,y)> dA


DISCRETISATION in 2D

The domain D has to be discretised in N nodes and T triangles. For linear elements the Nodes coincide with vertices, but for 
quadratic elements there will Nodes defined on the mid-point of edges of triangles, so there will be more Nodes than vertices.

We define a test function Phi_i per node N_i which will be used to linearly interpolate tangent vectors located at nodes to all 
triangles incident to the node N_i. 

Since we are in 2D, the tangent space at node N_i is the same as any other tangent space everywhere in the domain. Therefore 
parallel transport is just the identify function or “teleportation”. 

The local basis for N_i is defined as T_v = (e_1, e_2), which is an orthonormal frame. The local basis at triangle T_f = (e_1, e_2), 
the connection R_i^f = Identity. The vectors at U(x,y) are expressed in terms of “global basis” (e_1, e_2)

The connection that parallel transport vectors from node N_i to triangle f: Phi_i^f(X) is the identity function: 
Phi_i^f(X)  =  R_i^f X  =  X.  So It can be omitted.

We define scalar tests functions with compact support K_i(x, y) for every node N_i to be equal to 0 or 1 and linear everywhere else:

K_j(N_i) = 1   if j = i 
K_j(N_i) = 0   if i != j
K_j(x, y) = linear between 0 and 1

We use those K_j to define the basis function Psi_i^f (see https://math.stackexchange.com/questions/1320672/system-of-equations-for-vector-valued-functions-problems):

Psi_i^f(P) = K_j(P) X(P)

Where X(P) is vector located at point P expressed in basis (e_1, e_2)

Now we discretise the solution U(x,y) to be a finite sum:

U(x,y) = sum_j^N( U_j Psi_j(x,y) )

We discretise the test functions V(x, y) = Psi_j(x,y) for every node N.

In this case the test function is the same as the shape function Psi_j, when that happens this is the so called Galerkin FEM. 

Replacing first the test function V(x,y) with Psi_j(x,y) we get a system of integral equations, one equation per Node/test 
function Psi_j(x,y):

Integral_Di Nabla(Psi_i(x,y)) : Nabla(U(x,y) dx = Integral_D <F(x,y), Psi_j(x,y)> dA

where the integrals are now per triangle-set around node N_i or test function Psi_i

We discretise the Nabla(V) : Nabla(U) for each element Ti on the triangle-set around node N_i as contribution of nodes N_j:

Nabla U_i(x,y) = sum_i U_i Nabla Psi_i(x,y)
Nabla V_i(x,y) = Nabla Psi_i(x,y)

Plugging into: Nabla(V) : Nabla(U) 
 
Nabla(V) : Nabla(U)  = < Nabla Psi_i(x,y), sum_i U_i Nabla Psi_i(x,y)>
                     = sum_j < Nabla Psi_i(x,y), U_j Nabla Psi_j(x,y)> 
                	 = sum_j U_j < Nabla Psi_i(x,y), Nabla Psi_j(x,y)> 

< Nabla V, Nabla U>_i = sum_j U_j <Nabla Psi_j, Nabla Psi_i>    (one equation with many unknowns U_j for each test function i)

Integral_Di( < Nabla V_i, Nabla U_j> dA)  = sum_j U_j  Integral_Di( <Nabla Psi_j, Nabla Psi_i> dA)

So we get the following system of linear equations:

sum_j U_j  Integral_Di( <Nabla Psi_j, Nabla Psi_i> dA) = Integral_D i( F Psi_i dA)

Which can be written as:

M_ij = Integral_Di( <Nabla Psi_j, Nabla Psi_i> dA)
F_i   = Integral_Di ( F Psi_i dA)

M_ij U_i = F_i

Which can be assemble to form a global matrix system (we choose to have shape functions Psi_i = test functions so we 
have same number of equations than unknowns)

M U = F


INTEGRATION

We have to solve:

Integral_Di( <Nabla Psi_j, Nabla Psi_i> dA)

Given a vector valued function F(x,y):

F(x,y) = [ K(x,y)  X_1(x,y),   K(x,y)  X_2(x,y) ]

Where K(x,y) is a scalar valued function (scalar field) and constant vector X(x,y) (constant vector field) the Gradient/Jacobian is:

d F(x,y)=   d/dx K(x,y)  X_1(x,y),    d/dy K(x,y)  X_1(x,y)
                  d/dx K(x,y) X_2(x,y),    d/dy K(x,y) X_2(x,y)

The above can be also written as the following tensor product:

d F(x,y)  = grad K(x,y) @ X(x,y)

Nabla Psi_j(X(P)) = Nabla( K_j(P) X(P) )
                               = grad (K_j(P)) @ X(P) + K_j(P) Nabla(X(P))   // Leibnitz rule


Since X(P) is a constant vector field over the triangle then Nabla(X(P)) = 0

Nabla Psi_j(X(P))  = grad (K_j(P)) @ X(P)

Nabla Psi_j(X(P))  = grad (K_j(P)) @ (X_1 e_1, X_2 e_2)

Nabla Psi_j(X(P))  = [ grad (K_j(P)) @ (e_1, e_2) ] X


(grad_i  Kronecker  X) =    grad_1  * e_1     grad_1  * e_2       X_1         =  grad_1  * (e_1, e_2)  X
                                             grad_2 * e_1     grad_2 * e_2      X_2             grad_2 * (e_1, e_2)  X
2x1                           1x2       2x2                                                  2x1            2x1


Abstracting the X we get:

Nabla Psi_j  = grad (K_j(P)) @ (e_1, e_2)

So the inner product of G_i with G_j gives the following 2x2 matrix:


 G_i^T = (e_1, e_2) grad (K_i(P))^T
 G_j     = grad (K_j(P)) (e_1, e_2)^T 

 G_i^T G_j = (e_1, e_2) grad (K_i(P))^T grad (K_j(P)) (e_1, e_2)^T 

 G_i^T G_j = <grad_i, grad_j> (e_1, e_2) (e_1, e_2)^T 

(e_1, e_2) (e_1, e_2)^T  = Identity       // applying orthonormality


Seo we get:

G_i^T G_j = <grad_i, grad_j> Identity_2x2 

Later on we will use the above to define the Gaussian Quadrature of the shape function as:

w_ij = area_k * <grad_i(P), grad_j(P)> Identity_2x2

Where w_ij is symmetric I.e., w_ij == w_ji


DISCRETISATION in 3D

The domain D has to be discretised in N nodes and T triangles. For linear elements the Nodes coincide with vertices, but for 
quadratic elements there will Nodes defined on the mid-point of edges of triangles, so there will be more Nodes than vertices.

We define a test function Phi_i per node N_i which will be used to linearly interpolate tangent vectors located at nodes to 
all triangles incident to the node N_i. 

First we define a local basis of tangent T_v space per node N_i in which tangent vector of U(x,y) will reside. We will also 
need to define a local basis per triangle T_f which will be used to interpolate vectors inside triangles. Finally we need a 
“connection” R_i^f which parallel transport a vector in tangent space at node N_i basis to a triangle local basis:

R_i^f = T_f^T    Q_vf      T_v
      2x3        3x3        3x2
where:

R_i^f = is the connection at node N_i restricted to triangle f (2x2 matrix). Parallel transport a vector in tangent space 
at node N_i basis to a triangle local basis
T_f  = is the 3x2 local basis for triangle f : T_f = [e_i, e_j]^T
Q_vf = is the 3x3 rotation matrix aligning node N_i tangent plane with face f plane.
T_v = is the 3x2 local basis for tangent space at node N_i: T_v = [e_i, e_j]^T

We the above connection define the mapping Phi_i^f at node N_i restricted to triangle f:

Phi_i^f(X) = R_i^f X

Which parallel transport a vector X in tangent space basis at node N_i basis to a triangle local basis.

We also define scalar tests functions with compact support K_i(x, y) for every node N_i to be equal to 0 or 1 and linear everywhere else:

K_j(N_i) = 1   if j = i 
K_j(N_i) = 0   if i != j
K_j(x, y) = linear between 0 and 1

We use those K_j to define the basis function Psi_i^f:

Psi_i^f(P) = K_j(P) Phi_i^f(X) = K_j(P) R_i^f X

Where P is node’s position or any position inside triangle and X(P) is tangent vector at position P. Note that while P can vary 
inside a triangle, X(P) is constant inside a triangle.

Now we discretise the solution U(x,y) to be a finite sum:

U(x,y) = sum_j^N( U_j Psi_j(x,y) )

We discretise the test functions V(x, y) = Psi_j(x,y) for every node N.

In this case the test function is the same as the shape function Psi_j, when that happens this is the so called Galerkin FEM. 

Replacing first the test function V(x,y) with Psi_j(x,y) we get a system of integral equations, one equation per Node/test 
function Psi_j(x,y):

Integral_Di Nabla(Psi_i(x,y)) : Nabla(U(x,y) dx = Integral_D <F(x,y), Psi_j(x,y)> dA

where the integrals are now per triangle-set around node N_i or test function Psi_i

We discretise the Nabla(V) : Nabla(U) for each element Ti on the triangle-set around node N_i as contribution of nodes N_j:

Nabla U_i(x,y) = sum_i U_i Nabla Psi_i(x,y)
Nabla V_i(x,y) = Nabla Psi_i(x,y)

Plugging into: Nabla(V) : Nabla(U) 
 
Nabla(V) : Nabla(U)  = < Nabla Psi_i(x,y), sum_i U_i Nabla Psi_i(x,y)>
                                    = sum_j < Nabla Psi_i(x,y), U_j Nabla Psi_j(x,y)> 
                                    = sum_j U_j < Nabla Psi_i(x,y), Nabla Psi_j(x,y)> 

< Nabla V, Nabla U>_i = sum_j U_j <Nabla Psi_j, Nabla Psi_i>    (one equation with many unknowns U_j for each test function i)

Integral_Di( < Nabla V_i, Nabla U_j> dA)  = sum_j U_j  Integral_Di( <Nabla Psi_j, Nabla Psi_i> dA)

So we get the following system of linear equations:

sum_j U_j  Integral_Di( <Nabla Psi_j, Nabla Psi_i> dA) = Integral_D i( F Psi_i dA)

Which can be written as:

M_ij = Integral_Di( <Nabla Psi_j, Nabla Psi_i> dA)
F_i   = Integral_Di ( F Psi_i dA)

M_ij U_i = F_i

Which can be assemble to form a global matrix system (we choose to have shape functions Psi_i = test functions so we have 
same number of equations than unknowns)

M U = F

COVARIANT DERIVATIVE Nabla Psi_i 

We use the Cartan’s moving frames formalism. Moving frames are orthonormal frames attached to each point of surface. Those 
frames lie in the tangent space associated to each point, thus the frames are orthogonal unit tangent vectors. Those frames don’t form a holonomic basis i.e., those are not vector fields associated to a directional derivative or partial derivative operator. Moving frames are not continuous over the surface and its length is always 1 (holonomic basis can have different lengths, thus their dot product defines the metric tensor). For this reason the Crhistoffel Symbols is a concept associated to differentiation of holonomic basis while Connection Form is a concept associated with “differentiation” of moving frames. However in the flat space under cartesian coordinates the Christoffel Symbols and Connection Form are the same.

First recap definition of covariant derivative of a vector field U in direction of vector field w (https://en.wikipedia.org/wiki/Connection_form#Example:_the_Levi-Civita_connection):

E = (e_1, e_2)             // local frame field I.e., moving frame seen as a vector
U = U_1 e_1 + U_2 e_2      // vector field expressed in moving frames as sum_i U_i e_i 
w = w_1 e_1 + w_2 e_2      // vector field expressed in moving frames as sum_i w_i e_i  

W_ij(e_k)  =  W_11(e_k)        W_12(e_k)         // Matrix of 1-forms    a.k.a. connection form
                      W_21(e_k)       W_22(e_k)

Grad U_i = gradient of scalar component U_i which is a vector (a, b) = a e_1 + b e_2.  Where a = d/de_1 U_i   and   b= d/de_2 U_i

Nabla(U, w) = w_1 Nabla(U, e_1) + w_2 Nabla(U, e_2)       // Nabla(U, w) = <Nabla(U), w>

Nabla(U_i, e_j) = <Grad U_1, e_j> e_1  + ( W_11(e_j)  U_1 + W_12(e_j)  U_2 ) e_1
                              <Grad U_2, e_j>e_2 + ( W_21(e_j) U_1 + W_22(e_j) U_2 ) e_2
                                                                
Nabla(U, e_j) = <Grad U, e_j> E  + <W(e_j) U,  E>

Putting back w = w_1 e_1 + w_2 e_2

Nabla(U, w) = <Jacobian(U) w, E>  + <W(w) U,  E>

The above expression coincides with definition in Liu, Tong, De Goes, Desbrun “Discrete Connection and Covariant 
Derivative for Vector Field Analysis and Design”. 

<Jacobian(U) w, E>= The multiplication of Jacobian(U) with w is giving the directional derivative.
<W(w) U,  E> = The multiplication of W(w) with U is parallel transporting U to the tip of w

We can express Nabla(U, w) just in terms of U without passing w as argument.

Nabla(U_i) = (Grad U_i @ (e_1, e_2) + W_kji U_i) w_k

Nabla(U, w) = Jacobian(U) w  + W_ijk U_j w_k
Nabla(U, w) = [ Jacobian(U) + W(U) ]  w
Nabla(U) = Jacobian(U) + W(U)

With the above representation we can contract Nabla(U) with vector w to get the classic Levi-Civita connection:

Nabla(U, w) = <Nabla(U), w>

So we can rewrite Nabla(U) without taking into account w as (see https://math.stackexchange.com/questions/3468914/covariant-derivative-of-a-vector-field):

Nabla(U) = Jacobian(U) + T_ijk U_k

With that definition we define Nabla(Psi_i) as:

Nabla(Psi_i) = Nabla(K_j(P) Phi_i X(P))
Nabla(Psi_i) = Nabla(K_j(P)) Phi_i X(P) + K_j(P) Nabla(Phi_i (X(P)))     //Leibnitz rule

Nabla(Phi_i (X)) = Nabla(R_i^f X)
Nabla(R_i^f X) = Jacobian(R_i^f X) + T_ijk (R_i^f X)_k = 0

Nabla(R_i^f X) = 0 because:
Jacobian(R_i^f X) = 0 as R_i^f X is constant function inside triangle (constant vector field on flat triangle)
T_ijk (R_i^f X)_k = 0 as Connection Forms T_ijk are zero for orthonormal euclidean coordinates in flat space, and this 
is flat space since R_i^f X is a tangent vector in the face plane where vector fields are constant parallel vector fields.


Nabla Psi_i = Nabla K_j(P) Phi_i (X)
            = Grad K_j(P) @ Phi_i (X) + K_j(P) * Nabla_X Phi_i (X)


Nabla Psi_i = Grad K_j(P) @ Phi_i (X)       //since  Nabla_X Phi_i (X) = 0

where  
Phi_i (X) = is the vector at node (in node’s tangent space) mapped to face tangent space as 2x1 column
grad_X K_j(P) = is a 2x1 gradient of shape function K_j(P) in face basis as columns
K_j(P) = scalar-valued shape function with compact support


Nabla Psi_i = grad K_j(P) @ Phi_i (X)

To clarify the use of tensor “dyadic” product to express Jacobian matrix check the following:

K_j(X) V = | K_j(X) Vx,      K_j(X) Vy | 

d/X K_j(X) V = | d/x K_j(X) Vx,      d/x K_j(X) Vy |
                         | d/y K_j(X) Vx,      d/y K_j(X) Vy |
                     = | d/x K_j(X), d/y K_j(X) | @ | Vx, Vy | 
                     = Grad_X K(X) @ V

Expressing this in matrix notation we get:

Nabla Psi_i  = grad K_j(P) * Phi_i(X)^T
2x2                             2x1              1x2

It is convenient to abstract the argument X from above equation. For that we express it in tensor form as:

G_ij = grad_i @ R_jk X_k  //tensor product of two vectors grad_i with R_j

Abstracting X_k:

G_ijk = (grad_i @ R_jk)    // we get an rank-3 tensor

It is known that the contraction of a rank-3 tensor with a rank-1 tensor can be expressed in matrix form using the Kronecker product as:

(grad_i  Kronecker R_jk) =  grad_1  * R_ik     X_1         =  grad_1  * R X
                                               grad_2 * R_ik     X_2            grad_2 * R X
2x1                          2x2         4x2                     2x1               4x1

Where the resulting 4x1 vector (a b c d)^T is isomorphic to the 2x2 matrix:

 (a b c d)^T     isomorphic to      a b
                                                      c d

We will call G (covariant gradient) the 4x2 matrix:

G_i = (grad_i  Kronecker R_jk)
4x2

So the inner product of G_i with G_j gives the following 2x2 matrix:

 G_i^T G_j = <grad_i, grad_j> R_i^T R_j

We we can check:

grad_1 R^T  grad_2 R^T    grad_1 R    =   (grad_1*grad_1 + grad_2*grad_2) R^T R
                                             grad_2 R      


Later on we will use the above to define the Gaussian Quadrature of the shape function as:

w_ij = area_k * <grad_i(P), grad_j(P)> R_i^T R_j

Where w_ij is NOT symmetric, I.e., w_ij != w_ji

Definition of Shape Functions and Test Functions on Elements

Assembly

The process of creating the matrix M_ij is called assembly. It consist of creating the matrix M^e_ij for each element 
(triangle) and then add it to the global M_ij matrix. 

The global matrix M_ij is 2Nx2N where N is number of mesh nodes. In theory the matrix M^e_ij is a super-sparse version 
of M_ij just with a few non-zero entries which corresponds to a single element. So

M_ij = sum_e M^e_ij

In practice the M^e_ij is not super-sparse matrix of size 2Nx2N, but it is only a 6x6 matrix (in case of linear triangle elements)

But there is a map from indices from the 6x6 matrix to the 2Nx2N matrix as follows.

Each triangle is defined as having three nodes, let’s say triangle T_44 is made of nodes N_14, N_25 and N_31 of the mesh. 
Then we establish a local indexing of triangle with fixed nodes n_1, n_2 and n_3 and mapping 

n_1 = N_14
n_2 = N_25
n_3 = N_31

Then the 6x6 matrix would be like (keep in mind that Integral_Di( <Nabla Psi_14, Nabla Psi_14> dA) is a 2x2 matrix):

M_{0,0} = Integral_Di( <Nabla Psi_14, Nabla Psi_14> dA)       //2x2 block
M_{0,2} = Integral_Di( <Nabla Psi_14, Nabla Psi_25> dA)      //2x2 block
M_{0,4} = Integral_Di( <Nabla Psi_14, Nabla Psi_31> dA)       //2x2 block
…
M_{4,4} = Integral_Di( <Nabla Psi_31, Nabla Psi_31> dA)       //2x2 block

M^e_ij is not symmetric so we need to compute the full combinations.

Then the matrix needs to be mapped back to global numbers:

M_{2*14    , 2*14}     = M^e_(0, 0)
M_{2*14    , 2*14+1} = M^e_(0, 1)
M_{2*14+1, 2*14}     = M^e_(1, 0)
M_{2*14+1, 2*14+1} = M^e_(1, 1)

M_{2*14    , 2*25}     = M^e_(0, 2)
M_{2*14    , 2*25+1} = M^e_(0, 3)
M_{2*14+1, 2*25}     = M^e_(1, 2)
M_{2*14+1, 2*25+1} = M^e_(1, 3)

M_{2*14.   , 2*31}     = M^e_(0, 4)
M_{2*14    , 2*31+1} = M^e_(0, 5)
M_{2*14+1, 2*31}     = M^e_(1, 4)
M_{2*14+1, 2*31+1} = M^e_(1, 5)
…
M_{2*31    , 2*31}     = M^e_(4, 4)
M_{2*31    , 2*31+1} = M^e_(4, 5)
M_{2*31+1, 2*31}     = M^e_(5, 4)
M_{2*31+1, 2*31+1} = M^e_(5, 5)


Reference Elements

One of the novelties of FEM is to be able to compute the integrals per individual elements on an easy and generic way by 
using so called Reference Elements.

For triangle elements we just need a single reference triangle

P3
º
|    \
|          \
|                 \
º——————º
P1               P2


P1 = (0,0)
P2 = (1,0)
P3 = (0,1)

N_1(u,v) = 1 - u - v
N_2(u,v) = u
N_3(u,v) = v

| x |  = | x_2 - x_1   x_3 - x_1  |   | u | + | x_1 |
| y |     | y_2 - y_1    y_3 - y_1  |   | v |    | y_1 |
| z |     | z_2 - z_1   z_3 - z_1  |             | z_1 |
3x1      3x2                                  2x1    3x1


| x - x_1 |  = | x_2 - x_1   x_3 - x_1  |   | u |
| y - y_1 |     | y_2 - y_1    y_3 - y_1  |   | v |
| z - z_1 |     | z_2 - z_1   z_3 - z_1  |   
3x1              3x2                                   2x1

X = B U + X_1

(B^T B)^-1 B^T X = U
B^-* X = U

where B^-* = (B^T B)^-1 B^T

But 

F(U)       = B U + X_1                              In K
F^-1(X)  = (B^T B)^-1 B^T (X - X_1)   In Ref 
F^-1(X)  = B^-* (X - X_1)                      where B^-* = (B^T B)^-1 B^T

N^k_i(X) = N_i( F^-1(X) ) = N_i º F^-1(X)

grad_X N^k_i(X) = grad_X N_i( F^-1(X) )
                              = grad_U N_i( F^-1(X) ) * grad_X F^-1(X)    By chain rule
                              = B^-*T grad_U N_i( F^-1(X) )

grad_U N_i(U) = grad_U N^k_i( F(U) )
                         = grad_X N^k_i( F(U) ) * grad_U F(U)     By chain rule
                         = B^*T grad_X N^k_i( F(U) )

grad_U N_1(U) = (-1, -1)
grad_U N_2(U) = (1, 0)
grad_U N_3(U) = (0, 1)

grad_X F^-1(X) = grad_X B^-* (X - X_1) 
                           = grad_X B^-* X - grad B^-* X_1     (but grad B^-* X_1 = 0 as there is no X)
                           = grad_X B^-* X
                           = B^-*T      (transposed, where B^-*T is a 3x2 matrix)


grad_U N_1( F^-1(X) ) = (-1, -1)
grad_U N_2( F^-1(X) ) = (1, 0)
grad_U N_3( F^-1(X) ) = (0, 1)

grad_X N^k_1(X) = B^-T (-1, -1)
grad_X N^k_2(X) = B^-T (1, 0)
grad_X N^k_3(X) = B^-T (0, 1)

w_ij = integral_k <Nabla Psi_i(X), Nabla Psi_j(X)> dA

Since:

Nabla Psi_i(P) = Grad N^k_i(P) tensor Phi_i(X(P))

We get:

w_ij = integral_k <Grad N^k_i(P) tensor Phi_i(X(P)), Grad N^k_j(P) tensor Phi_j(X(P))> dA

I^k_ij(X) = <Grad N^k_i(P) tensor Phi_i(X(P)), Grad N^k_j(P) tensor Phi_j(X(P))>

w_ij = area K / 3 ( I^k_ij(X_1) +  I^k_ij(X_2) +  I^k_ij(X_3))

Which is an exact quadrature formula for linear function on a triangle. Since grad N^k_i(P) does not depend on P I.e., it is 
constant, and also Phi_i(X(P)) basis function is constant per triangle the quadrature expression reduces to:

w_ij = area K <Grad N^k_i(P) tensor Phi_i(X(P)), Grad N^k_j(P) tensor Phi_j(X(P))>

We can define this as:

w_ij = area_k * <grad_i(P), grad_j(P)> R_i^T R_j

Where R_i and R_j are the connection matrices of nodes N_i and N_j respectively which “parallel transport” tangent vectors from 
nodes to the triangle and allows to work with a “flat” constant metric inside triangle induced by R^3 metric.
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

