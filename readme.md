# FreeAlgebra
## An autodifferentiable libary for symbolic algebra

This project centers on defining and manipulating **finitely generated algebras** in a differentiable manner, particularly within the context of quantum and mathematical physics. Key features:
- **Custom algebra definition**: Define algebras via commutation relations, specify traces; `FreeAlgebra` can infer the algebra's multiplication rules automatically.
- **Automorphism-based algebra construction**: Construct new finitely-generated algebras by specifying automorphisms, with automatic inference of new commutation relations and multiplication rules. 
- **Tensor products and expansion**: Expand algebraic systems by constructing tensor products or finite powers, useful for modeling multi-particle or iterative systems.
- **Differentiability**: all operations are built on top of Pytorch with support for automatic differentiation. 

This codebase is well-suited for advanced research in quantum information, algebraic geometry, and the computational modeling of algebraic structures in physics. It was originally developed to analyze Grassmann integration and the Clifford-Grassmann Fourier transform, and is not being actively maintained. 

### Files and Descriptions

#### 1. `complex.h`
- **Complex Numbers and Differentiability**: Implements complex numbers, including arithmetic and conjugation operations. It integrates with PyTorch to provide differentiable computations, enabling optimization and learning in quantum contexts.

#### 2. `Elements.h`
- **Algebraic Element Representation**: Defines the `AlgebraElement` class for representing elements of an algebra.
- **Core Operations**: Implements basic arithmetic (addition, multiplication) and algebraic properties (conjugation, scalar multiplication).
- **Validation and Filtering**: Ensures elements maintain consistency with the algebra’s defined relations.

#### 3. `BaseAlgebraRelations.h`
- **Commutation Relations**: Allows custom algebras to be defined by specifying commutation relations between generators.
- **Multiplication Rules**: The multiplication rules are inferred based on these commutation relations.
- **Trace Specification**: Optionally defines the trace operation for algebraic elements, which is crucial for analyzing the algebra’s properties.

#### 4. `ProductPowerAlgebra.h`
- **Tensor Product of Algebras**: Defines tensor products of finitely generated algebras. This is used to model multi-particle quantum states or higher-dimensional algebraic operations.
- **Expansion of Algebras**: Allows construction of more complex algebras by combining simpler ones through tensor product operations.

#### 5. `Automorphism.h`
- **Automorphisms for Defining New Algebras**: Enables the specification of a new algebra by providing an automorphism of an existing algebra.
- **Inference of Commutation Relations**: New commutation relations are inferred based on the automorphism and the original algebra's structure, allowing efficient definition of related algebras.

#### 6. `Fermions.h`
- **Fermionic Systems**: Implements algebraic representations specific to fermions, focusing on anti-commutation relations.
- **Fermionic Algebra Constructs**: Provides abstractions for defining and manipulating fermionic algebras, which are critical in quantum field theory and many-body quantum physics.

#### 7. `FinitePowerAlgebras.h`
- **Finite Power Algebras**: Implements functionality for constructing finite power algebras, which involve repeated multiplication of algebraic elements.
- **Iterative Constructs**: Useful for iterating algebraic structures, such as powers of generators or interaction terms in algebraic expansions.

#### 8. `utils.h`
- **Utility Functions**: Provides auxiliary functions for matrix operations, enumerations, and algebraic validation.
  - **Matrix Conversions**: Supports conversion of algebraic elements into matrices, particularly for compatibility with numerical libraries like Eigen.
  - **Enumeration and Validation**: Implements enumeration of algebraic terms and validation of element consistency, supporting algebraic structure exploration.