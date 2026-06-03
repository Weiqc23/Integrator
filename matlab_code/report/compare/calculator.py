from math import log10
import cmath
import math
x = 200*0.1445*log10(7/(0.011259))
b = 7.58/(log10(7/(13.9*1e-3)))*1e-6

print('x')
print(x)
print('b')
print(b)

S2_1 = (160+120j) + (160**2 + 120**2)/(220**2)*(1.12+28.9j)
print('S2_1')
print(S2_1)

S2_2 = 220*220* (6.612+54.63j)*1e-6
print('S2_2')
print(S2_2)

S2 = S2_1 + S2_2
print('S2')
print(S2)

S1 = S2 + abs(S2)**2 / 220**2 * (15.75 + 80.73j)
print('S1')
print(S1)

U2 = 240 - S1.conjugate() / 240 * (15.75 + 80.73j)
print('U2')
print(U2)

U3 = U2 - (S2_1/U2).conjugate() * (1.12 + 28.9j)
print('U3')
print(U3)


u3 = U3 * 121 / 220
print('u3')
print(u3)

eta = 161.2 / 176.7 * 100
print('eta')
print(eta)


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

u1 = 110 - (SA1/110).conjugate()*(3+18j)
print('u1')
print(u1)

u2 = 110 - (SA2/110).conjugate()*(4+15j)
print('u2')
print(u2)




Deq = 