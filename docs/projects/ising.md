# Transverse Field Ising Model

Numerical study of the transverse field Ising model (TFIM) for *n* spin-½ particles.  
The Hamiltonian is built and diagonalized using NumPy.

---

## Report
<iframe src="report.pdf" width="100%" height="700px"></iframe>
---

## Repository
[GitHub Repository](https://github.com/Tomas-Gradowski/ising-model)

---
## Code
??? note "Show full code"
    ```python
    # ising.py
    import numpy as np
    
    Sx = np.array([[0, 1],[1, 0]], dtype = complex)
    Sz = np.array([[1, 0],[0, -1]], dtype = complex)
    Id = np.array([[1, 0],[0, 1]], dtype = complex)

    def kron_list(lst):
        a = lst[0]
        for i in lst[1:]:
        a = np.kron(a, i)
      return a
    
    def sx(k, n):
        a = [Id]*n
        a[k] = Sx
      return kron_list(a)
    def sz(k, n):
        b = [Id]*n
        b[k] = Sz
    return kron_list(b)
    def ising(n):
      H = np.zeros((2**n, 2**n), dtype=complex)
      for k in range(n-1):
        H += sz(k, n) @ sz(k+1, n)
      for k in range(n):
        H -= sx(k, n)
    return H

    def printH(n):
      H = ising(n)
      print("The hamiltonian is:")
      print(np.real_if_close(H))
    def gndstate(n):
      H = ising(n)
      E0 = np.min((np.linalg.eigh(H))[0])
      print("its grounstate energy is:")
      print(E0)

    def memory(n):
      print("the required space in terabytes is:")
      a = (2**(2*n)*(np.zeros((1,1), dtype=complex).nbytes))/10**12
    print(a)
    ```
---
## Results
- For n = 2 → E₀ ≈ −2.236  
- For n = 10 → E₀ ≈ −12.38  
Energy per spin converges to ≈ −1.238.  
Memory grows as \(2^{2n}\), making diagonalization impractical for \(n \gtrsim 19\).

---

## Files
- `ising.py`: main script for Hamiltonian construction and diagonalization.  
- `report.pdf`: derivation and results summary.


