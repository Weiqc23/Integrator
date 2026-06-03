from math import log10
import cmath
import math

z_all = 3-18j + 4-8j + 4-15j
print('z_all')
print(z_all)



SA1 = (36+16j)*(8-23j)/z_all + (20+8j)*(4-15j)/z_all
print('SA1')
print(SA1)

SA2 = (36+16j)*(3-18j)/z_all + (20+8j)*(7-26j)/z_all
print('SA2')
print(SA2)

S12 = SA1 - (36+16j)
print('S12')
print(S12)

SA = (SA1) + abs(SA1)**2 / 110**2 * (3+18j)
print('SA')
print(SA)

U1 = 110 - (SA/110).conjugate() * (3+18j)
print('U1')
print(U1)

S2_1 = -S12 + abs(S12)**2 / 110**2 * (4+8j)
print('S2_1')
print(S2_1)

S2_2 = S2_1 + (20+8j)
print('S2_2')
print(S2_2)

SA_1 = S2_2 + abs(S2_2)**2 / 110**2 * (4+15j)
print('SA_1')
print(SA_1)

U2 = 110 - (SA_1/110).conjugate() * (4+15j)
print('U2')
print(U2)


U1_1 = U2 - (S2_1/U2).conjugate() * (4+8j)
print('U1_1')
print(U1_1)


k = (16-13.57) * (11**2 + 41**2)/ 110**2 / 41
print('k')
print(k)