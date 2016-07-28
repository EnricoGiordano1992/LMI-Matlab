mu = 0.01
Ms = 0.61
Bs = 11
Kv = 87.8
Kp = 4000

s = tf('s')

num = mu * Kp
den = Ms*s^2 + (Bs + Kv)*s + Kp

H = num / den

discrete_sys = c2d(H, 0.001)

ss(discrete_sys)
