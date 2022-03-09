# Hsiao-Hsuan(Charlie) Chin
# March 4, 2022
import math
import matplotlib.pyplot as plt
import cmath
from Parameters import *

# right hand side of phi equation
def rhs_dphi_dt(xx):
    return xx ** 1.5 / MM

# right hand side of xx equation
def rhs_dxx_dt(xx):
    t1 = chi1 + chi2
    t2 = t1 * t1 / 4.0
    t4 = t2 * t1 / 2.0
    t9 = 0.3141592653589793 * 10 * eta
    t11 = eta * eta
    t14 = chi1 - chi2
    t15 = delta * t14 / 2.0
    t17 = t14 * t14 / 4.0
    t19 = t17 * t14 / 2.0
    t20 = delta * t19
    t22 = t1 * eta / 2.0
    t24 = t1 * t11 / 2.0
    t26 = t17 * t1 / 2.0
    t28 = eta * t4
    t30 = t11 * eta
    t39 = delta * delta
    t40 = t39 * t17
    t43 = eta * t17
    t44 = t43 * t1 / 2.0
    t46 = t11 * t17
    t51 = -505.0 / 8.0 * t4 - 2529407.0 / 54432.0 * chi2 - 4415.0 / 4032.0 * 0.3141592653589793 * 10 \
          - 2529407.0 / 54432.0 * chi1 + 358675.0 / 6048.0 * t9 + 91495.0 / 1512.0 * \
          t11 * 0.3141592653589793 * 10 - 2529407.0 / 27216.0 * t15 - 505.0 / 8.0 * t20 + 10772921.0 / 54432.0 \
          * t22 - 398017.0 / 2016.0 * t24 - 523.0 / 8.0 * t26 + 979.0 / 24.0 * t28 + 2045.0 / 432.0 * t30 * t1 \
          - 11.0 / 24.0 * t11 * t4 + 12.0 * t17 * 0.3141592653589793 * 10 + 12.0 * t2 * 0.3141592653589793 * 10 \
          - 62.0 * t40 * t1 + 7007.0 / 24.0 * t44 - 3397.0 / 48.0 * t46 * t1 - 48.0 * t43 * 0.3141592653589793 * 10
    t52 = MM * MM
    t53 = 1 / t52
    t55 = m1 * m1
    t58 = t15 * eta
    t60 = delta * eta
    t61 = t60 * t19
    t63 = t11 * delta
    t68 = t15 * t2
    t71 = m2 * m2
    t74 = 1 / MM
    t75 = delta * t74
    t92 = t53 * eta
    t93 = chi1 * t55
    t96 = chi2 * t71
    t99 = t53 * t11
    t110 = t14 * t2 / 2.0
    t111 = t60 * t110
    t115 = 0.3141592653589793 * 10 * t1 / 2.0
    t118 = -1113425.0 / 13608.0 * t53 * chi1 * t55 + 845827.0 / 6048.0 * t58 + 742.0 / 3.0 * t61 \
           -41551.0 / 1728.0 * t63 * t14 + t63 * t19 / 8.0 - 1515.0 / 8.0 * t68 - 1113425.0 / 13608.0 * t53 * chi2 * \
           t71 - 11981.0 / 63.0 * t75 * eta * m1 * chi1 + 11981.0 / 63.0 * t75 * eta * m2 * chi2 + 3407.0 / 72.0 * \
           t75 * t11 * m1 * chi1 - 3407.0 / 72.0 * t75 * t11 * m2 * chi2 + 673643.0 / 1512.0 * t92 * t93 + \
           673643.0 / 1512.0 *t92 * t96 - 24829.0 / 216.0 * t99 * t93 - 24829.0 / 216.0 * t99 * t96 + 671.0 / 336.0 \
           * t75 * m1 * chi1 -671.0 / 336.0 * t75 * m2 * chi2 + 917.0 / 12.0 * t111 + 3.0 / 8.0 * t63 * t110 + \
           24.0 * t15 * t115
    v = math.sqrt(xx)
    t120 = v * v
    t121 = t120 * v
    t122 = t120 * t120
    t126 = 0.3141592653589793 * 10 * 0.3141592653589793 * 10
    t128 = math.log(v)
    t130 = math.log(2.0)
    t141 = -1712.0 / 105.0 * 0.5772156649015329 + 16.0 / 3.0 * t126 - 1712.0 / 105.0 * t128 \
           -3424.0 / 105.0 * t130 + 215.0 / 224.0 * t17 - 5605.0 / 2592.0 * t30 + 128495.0 / 2016.0 * t2 + 541.0 / \
           896.0 * t11 - 56198689.0 / 217728.0 * eta + 451.0 / 48.0 * eta * t126 - 2435.0 / 224.0 * t43 + 565.0 / \
           9.0 * t40
    t143 = eta * t2
    t148 = t15 * t1 / 2.0
    t150 = 0.3141592653589793E1 * delta
    t163 = 0.3141592653589793E1 * t53
    t171 = 16447322263.0 / 139708800.0 + 89.0 / 3.0 * t46 - 23441.0 / 288.0 * t143 + 1517.0 / \
           72.0 * t11 * t2 - 80.0 / 3.0 * t115 + 128495.0 / 1008.0 * t148 - 40.0 / 3.0 * t150 * t14 + 20.0 / 3.0 * \
           t9 * t1 + 31.0 / 6.0 * t150 * t74 * m1 * chi1 - 31.0 / 6.0 * t150 * t74 * m2 * chi2 - 16.0 * t163 * \
           t93 - 16.0 * t163 * t96 - 12733.0 / 576.0 * t60 * t14 * t1
    t191 = -31571.0 / 2016.0 * chi1 - 31571.0 / 2016.0 * chi2 - 3.0 / 4.0 * t4 + 9.0 / 4.0 * t111 \
           -31571.0 / 1008.0 * t15 - 3.0 / 4.0 * t20 + 5791.0 / 63.0 * t22 - 79.0 / 3.0 * t24 - 9.0 / 4.0 * t26 + \
           9.0 / 4.0 * t28 - 189.0 / 8.0 * t9 + 1165.0 / 24.0 * t58 + 3.0 / 4.0 * t61 - 9.0 / 4.0 * t68 + 27.0 \
           / 4.0 * t44 -4159.0 / 672.0 * 0.3141592653589793 * 10
    t215 = t122 * t122
    dv_dt = 32.0 / 5.0 * (1.0 + (t51 + t118) * t122 * t121 + (t141 + t171) * t122 * t120 + t191 * t122 * v + \
            (13661.0 / 2016.0 * eta + 59.0 / 18.0 * t11 + 34103.0 / 18144.0 + 81.0 / 16.0 * t17 + 81.0 / 16.0 * \
            t2 - 20.0 * t43 - t143 / 4.0 + 81.0 / 8.0 * t148) * t122 + (-113.0 / 24.0 * chi1 - 113.0 / 24.0 * \
            chi2 - 113.0 / 12.0 * t15 + 19.0 / 3.0 * t22 + 4.0 * 0.3141592653589793E1) * t121 + (-11.0 / 4.0 * eta \
            - 743.0 / 336.0) * t120) * eta * t215 * v * t74

    return 2 * math.sqrt(xx) * dv_dt


# right hand side of m1, decreasing function of hawking radiation
def rhs_dm1_dt(m):
    return -10**-6


def rhs_dm2_dt(m):
    return -10**-6


# going one time step while updating values in both phi and xx lists
def one_rk4(index, dt):
    k1phi = rhs_dphi_dt(xx[index])
    k1xx = rhs_dxx_dt(xx[index])
    k2phi = rhs_dphi_dt(xx[index] + 0.5 * dt * k1xx)
    k2xx = rhs_dxx_dt(xx[index] + 0.5 * dt * k1xx)
    k3phi = rhs_dphi_dt(xx[index] + 0.5 * dt * k2xx)
    k3xx = rhs_dxx_dt(xx[index] + 0.5 * dt * k2xx)
    k4phi = rhs_dphi_dt(xx[index] + dt * k3xx)
    k4xx = rhs_dxx_dt(xx[index] + dt * k3xx)

    phi[index + 1] = phi[index] + dt * (k1phi + 2 * k2phi + 2 * k3phi + k4phi) / 6
    xx[index + 1] = xx[index] + dt * (k1xx + 2 * k2xx + 2 * k3xx + k4xx) / 6

# waveform 22(quantum numbers)
def hhat22(xx, phi):
    t1 = 1 + xx * (-107 / 42 + 55 / 42 * eta)
    t2 = xx ** 1.5 * (2 * math.pi)
    t3 = xx ** 2 * (-2173 / 1512 - 1069 / 216 * eta + 2047 / 1512 * eta ** 2)
    t4 = xx ** 2.5 * ((-107 / 21 + 34 / 21 * eta) * math.pi - 24 * complex(0,1) * eta)
    t6 = -856 / 105 * gamma - 1712 / 105 * math.log(2) - 428 / 105 * math.log(16 * xx)
    t7 = -(278185 / 33264 - 41 / 96 * math.pi ** 2) * eta - 20261 / 2772 * eta ** 2
    t8 = 114635 / 99792 * eta ** 3 + 428 * complex(0, 1) * math.pi / 105
    t5 = xx ** 3 * (27027409 / 646800 + 2 / 3 * math.pi ** 2 + t6 + t7 + t8)
    t9 = xx ** 3.5 * (-2173 * math.pi / 756 + (-2495 * math.pi / 378 + 14333 * complex(0, 1) / 162) * eta + \
                      (40 * math.pi / 27 - 4066 * complex(0, 1) / 945) * eta ** 2)
    orbit = t1 + t2 + t3 + t4 + t5 + t9


    c1 = xx ** (0.5) * (-2 * eta) * (chia ** 2 - chis ** 2)
    c2 = xx ** (1.5) * (-4 * delta * chia / 3 + 4 / 3 * (eta - 1) * chis + c1)
    c3 = xx ** (2.0) / MM ** 2 *(MM * (m1 * chi1 + m2 * chi2)) ** 2 / MM ** 2
    npspin = c2 + c3

    p1 = -0.5 * iota ** 2
    p2 = xx ** (0.5) * iota * delta * (1/3) * cmath.exp(complex(0, 1) * phi)
    p3 = delta * (chis - complex(0, 1) * chis)
    p4 = xx * (-0.5) * cmath.exp(complex(0, 1) * (alpha + phi)) * (chia - complex(0, 1) * chia + p3)
    ppspin = p1 + p2 + p4

    hhat = orbit + npspin + ppspin
    return hhat

# waveform 21
def hhat21(xx, phi):
    o1 = xx ** 0.5 +xx ** 1.5 * (-17 / 28 + 5 * eta / 7)
    o2 = xx ** 2.0 * (math.pi - complex(0, 1) / 2 - 2 * complex(0, 1) * math.log(2))
    o3 = xx ** 2.5 * (-43 / 126 - 509 * eta / 126 + 79 * eta ** 2 / 168)
    o5 = eta * (-995 / 84 - 3 * math.log(2) / 7) + 17 * math.log(2) / 14
    o4 = xx ** 3.0 * (-17 * math.pi / 28 + 3 * math.pi * eta / 14 + complex(0, 1) * (17 / 56 + o5))
    orbit = o1 + o2 + o3 + o4

    orbit = orbit * complex(0, 1) / 3 * delta

    n1 = complex(0, -1) * xx *(0.5) * (chia+ delta * chis)
    n3 = 1 / 42 * (-79 + 139 * eta) * MM * (m2 * chi2 - m1 * chi1)
    n2 = complex(0, -1) * xx ** 2 / MM ** 2 * (-43 / 21 * delta * (m1 ** 2 * chi1 + m2 ** 2 * chi2) + n3)
    npspin = n1 + n2

    p2 = 0.25 * iota ** 3 * (-5 / 3 * cmath.exp(complex(0, -1) * phi) - cmath.exp(3 * complex(0, -1) * phi))
    p1 = (iota * cmath.exp(complex(0, -1) * phi) + p2)
    p3 = xx ** 0.5 * delta * iota ** 2 * (5 / 12 - 0.25 * cmath.exp(2 * complex(0, -1) * phi))
    p6 = cmath.exp(complex(0, -1) * alpha) * (-1 + cmath.exp(2 * complex(0, -1) * phi))
    p5 = 0.25 * (chia - complex(0, -1) * chia + delta * (chis - complex(0, -1) * chis)) * p6
    p4 = xx * iota *(cmath.exp(complex(0, -1) * phi) / 42 * (-107 + 55 * eta) + p5)
    p8 = (1 - eta / 2) * (chis - complex(0, -1) * chis)
    p10 = delta *(chia + complex(0, -1) * chia) + (1 + 5 * eta / 6) * (chis + complex(0, -1) * chis)
    p9 = cmath.exp(complex(0, -1) * (alpha + phi)) * p10
    p7 = xx ** 1.5 *(-cmath.exp(complex(0, -1) * (alpha + phi)) * (delta * (chia - complex(0, -1) * chia) + p8) + p9)

    ppspin = complex(0, -1) * (p1 + p3 + p4 + p7)
    hhat = orbit + npspin + ppspin
    return hhat

# waveform 20
def hhat20(xx, phi):
    orbit = (-5 / (14 * math.sqrt(6)))
    npspin = 0

    p1 = 2 * iota ** 2 * math.cos(2 * phi)
    p2 = xx ** 0.5 * 4 * complex(0, 1) / 3 * iota * delta * math.sin(phi)
    p6 = (chis - complex(0, 1) * chis)
    p3 = xx / 3 * -cmath.exp(complex(0, 1) * (alpha + phi)) * (chia - complex(0, 1) * chia + delta * p6)
    p7 = (chia + complex(0, 1) * chia + delta * (chis+complex(0,1) * chis))
    p4 = cmath.exp(complex(0,-1) * (alpha + phi)) * p7
    p5 = xx * (-4 * complex(0, 1) / 3 * iota * math.sin(phi) * (chia + delta * chis))

    ppspin = p1 + p2 + p3 + p4 + p5
    ppspin = ppspin * 0.5 * math.sqrt(3 / 2)

    hhat = orbit + npspin + ppspin
    return hhat

# waveform 30
def hhat30(xx, phi):
    orbit = -2 / 5 * complex(0, 1) * math.sqrt(6 / 7) * xx ** 2.5 * eta
    npspin = 0
    ppspin = xx ** 0.5 * (-1 / (2 * math.sqrt(42))) * delta * iota * math.cos(phi)
    hhat = orbit + npspin + ppspin
    return hhat

# waveform 31
def hhat31(xx, phi):
    o1 = xx ** 0.5 + xx ** 1.5 * (-8 / 3 - 2 * eta / 3)
    o2 = xx ** 2.0 * (math.pi - 7 * complex(0, 1) / 5 - 2 * complex(0, 1) * math.log(2))
    o3 = xx ** 2.5 *(607 / 198 - 136 * eta / 99 + 247 * eta ** 2 / 198)
    o5 = eta * (-1 / 15 + 7 * math.log(2) / 3)
    o4 = xx ** 3.0 *(-8 * math.pi / 3 - 7 * math.pi * eta / 6 + complex(0, 1) * (56 / 15 + 16 * math.log(2)/3 + o5))
    orbit = o1 + o2 + o3 + o4
    orbit = orbit * complex(0, 1) * delta / (12 * math.sqrt(14))

    npspin = 0

    p2 = 3 * cmath.exp(2 * complex(0, 1) * phi)
    p1 = xx ** 0.5 * iota ** 2 * delta / 2 * (-11 / 2 + 135 / 2 * cmath.exp(-2 * complex(0, 1) * phi) - p2)
    p3 = xx * iota * 20 * (-1 + 3 * eta) * cmath.exp(complex(0, -1) * phi)
    p4 = xx ** 1.5 *(-16 * eta) * cmath.exp(complex(0, -1) * (alpha + phi)) * (chis + complex(0, 1) * chis)
    ppspin = p1 + p3 + p4
    ppspin = ppspin * complex(0, -1) * (-1 / (12 * math.sqrt(14)))

    hhat = orbit + npspin + ppspin
    return hhat



def main():
    global m1
    global m2
    global M
    global mu
    global eta
    global delta

    hr = []
    hi = []
    for i in range(len(t) - 1):
        m1 = m1 + rhs_dm1_dt(m1) * dt
        m2 = m2 + rhs_dm2_dt(m2) * dt

        M = m1 + m2
        mu = m1 * m2 / MM
        eta = mu / MM
        delta = (m1 - m2) / (m1 + m2)

        one_rk4(i, dt)
    for i in range(len(t)):
        z = hhat22(xx[i], phi[i])
        hr.append(z.real)
        hi.append(z.imag)

        # column 1: time  column 2: phi(t) column 3: velocity(t) column 4: real part of waveform
        print(t[i], phi[i], xx[i], hr[i])

    plt.plot(t,hr,"red")
    plt.plot(t,hi)
    plt.show()

main()