#Code which investigates and implements the qiskit HHL algorithm functionality for a 2x2 system, and 
#uses the observables feature to check the output against a classical linear solver.

import numpy as np
from qiskit.algorithms.linear_solvers.numpy_linear_solver import NumPyLinearSolver
from qiskit.algorithms.linear_solvers.hhl import HHL
from qiskit.algorithms.linear_solvers.matrices.tridiagonal_toeplitz import TridiagonalToeplitz
from qiskit.quantum_info import Statevector
from qiskit.algorithms.linear_solvers.observables import AbsoluteAverage, MatrixFunctional
from scipy.sparse import diags

#First two algorithms take numpy arrays as arguments (loses quantum advantage as exp(iAt) must be prepared)
matrix = np.array([[1, 1/2], [1/2, 1]]) #this is matrix for my problem
vector = np.array([0, 1]) #vector for my problem
naive_hhl_solution = HHL().solve(matrix, vector) #HHL 
classical_solution = NumPyLinearSolver().solve(matrix, vector / np.linalg.norm(vector)) #Classical solver (normalised b vector)
#We have a toeplitz matrix, so can pass as quantum circuit not numpy array
tridi_matrix = TridiagonalToeplitz(1, 1, 1 / 2)
tridi_solution = HHL().solve(tridi_matrix, vector)
#Look at classical solution, indeed the correct answer x
print('classical state:', classical_solution.state)
#The other two algorithms only give the state of the quantum system, circuit which generates this accessed below
print('naive state:')
print(naive_hhl_solution.state)
print('tridiagonal state:')
print(tridi_solution.state)
#Look at euclidean norm of the three outputs, all very close
print('classical Euclidean norm:', classical_solution.euclidean_norm)
print('naive Euclidean norm:', naive_hhl_solution.euclidean_norm)
print('tridiagonal Euclidean norm:', tridi_solution.euclidean_norm)
#Want to check that solution vectors from quantum algorithms are the same as the classical
naive_sv = Statevector(naive_hhl_solution.state).data
tridi_sv = Statevector(tridi_solution.state).data
# Extract the right vector components. 1000 corresponds to the index 8 and 1001 corresponds to the index 9
naive_full_vector = np.array([naive_sv[8], naive_sv[9]])
tridi_full_vector = np.array([tridi_sv[8], tridi_sv[9]])
naive_full_vector = np.real(naive_full_vector) #Remove complex parts added by computer inaccuracy
tridi_full_vector = np.real(tridi_full_vector)
print('naive raw solution vector:', naive_full_vector)
print('tridi raw solution vector:', tridi_full_vector)
#Divide by norms and multiply by euclidean norm
print('full naive solution vector:', naive_hhl_solution.euclidean_norm*naive_full_vector/np.linalg.norm(naive_full_vector))
print('full tridi solution vector:', tridi_solution.euclidean_norm*tridi_full_vector/np.linalg.norm(tridi_full_vector))
print('classical state:', classical_solution.state)
#Naive solution vector and classical solution both give the same answer! 


#Use observables feature, gives same answer for both quantum and classical - more useful for my problem
num_qubits = 1
matrix_size = 2 ** num_qubits
# entries of the tridiagonal Toeplitz symmetric matrix
a = 1
b = 1/2

matrix = diags([b, a, b], [-1, 0, 1], shape=(matrix_size, matrix_size)).toarray()
vector = np.array([0]*(matrix_size - 1)+[1])
tridi_matrix = TridiagonalToeplitz(1, a, b)

average_solution = HHL().solve(tridi_matrix, vector, AbsoluteAverage())
classical_average = NumPyLinearSolver().solve(matrix, vector / np.linalg.norm(vector), AbsoluteAverage())

#First observable is absolute average, same in both cases
print('quantum average:', average_solution.observable)
print('classical average:', classical_average.observable)

#Second observable is x^{T}Bx where B is a toeplitz matrix (a,b;b,a), defined through MatrixFunctional(a,b)
observable = MatrixFunctional(1, 1 / 2) #defined matrix A from problem

functional_solution = HHL().solve(tridi_matrix, vector, observable)
classical_functional = NumPyLinearSolver().solve(matrix, vector / np.linalg.norm(vector), observable)

print('quantum functional:', functional_solution.observable) #gives answer of 4/3 in both cases
print('classical functional:', classical_functional.observable)