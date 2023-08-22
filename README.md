# Compile Time LDL

Library for $LDL^T$ factorization of symmetric matrices whose sparsity is known at compile-time.
Intended for very small matrices or matrices with very small repeating blocks.
Avoids all looping and indexing overheads completely by unrolling all iterations of the factorization and generating very explicit computation steps.
This can improve runtime performance significantly but comes with much longer compilation times and possibly larger binary size.

The library is in early development and the API is subject to change.
