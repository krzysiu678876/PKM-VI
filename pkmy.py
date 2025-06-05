# Założenia projektowe
n_proj = 1200  # obr/min
F = 250  # N
Q = 100  # N
c = 260  # mm
d = 300  # mm
L10h = 110000  # obrotów
P_proj = 450  # Watt
epsilon = 0.02  # współczynnik do poślizgu

# Parametry Silnika
n_sil = 1400  # obr/min
P_sil = 0.55  # kW
M_sil = 3.750  # Nm

# Obliczenia
i = n_proj / n_sil * (1 - epsilon)  # przełożenie
print("przełożenie i =", i)
i_teo = n_proj / n_sil  # teoretyczne przełożenie
deltaI = (i - i_teo) / i * 100  # błąd przełożenia
print("błąd przełożenia i =", deltaI)

# Teraz mam średnicę koła podziałowego
D_1_min = 1200 * (P_sil / n_sil) ** (1 / 3)
print("D_1_min =", D_1_min)

# Daję średnicę koła napędowego jako tę z tabeli
D_1 = 90  # mm
D_2 = D_1 / i  # mm
print("średnica koła napędowego D_1 =", D_1)
print("średnica koła napędowego D_2 =", D_2)

# Teraz liczę siłę obwodową na D1 i D2
F_1 = 2 * 10**3 * M_sil / D_1
F_2 = 2 * 10**3 * M_sil / D_2

print("F_1 =", F_1)
print("F_2 =", F_2)
v_1 = 3.14 * D_1 * n_sil / (60 * 1000)  # m/s
print("v_1[m/s] =", v_1)

# Teraz liczę moce przenoszone na pasach
P_1 = F_1 * v_1
P_2 = F_2 * v_1
print("P_1 =", P_1)
print("P_2 =", P_2)

# Tutaj siła naciągu pasa
phi = 0.5  # kąt naciągu pasa
F_0 = F_1 / phi
print("F_0 =", F_0)

# Liczenie długości pasa
Ls = 980  # mm

# Liczenie siły naciągu pasa
# C25-materiał wału
Re = 280  # MPa
Rm = 460  # MPa
Zg = 200  # MPa
Zso = 130  # MPa
Zsj = 250  # MPa

# Liczenie odległości osi
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

aprim = np.sqrt(366**2 + 30**2)
print("odległość osi =", aprim)
alfa = 180 - 57 * (D_2 - D_1) / aprim
print("alfa =", alfa)

# Parametry wału
a = 0.33  # długość wału [m]
x2 = 0.0865  # pozycja pierwszej podpory [m]
x3 = 0.2335  # pozycja drugiej podpory [m]
d = 0.3  # pozycja siły Q [m]
x1 = 0.0225  # pozycja siły od pasa F0 [m]
Q = 100  # siła pionowa [N]
alfaprim = 4.47  # stopnia
Fx = F_0 * np.cos(alfaprim * np.pi / 180)  # siła naciągu pasa [N]
Fy = F_0 * np.sin(alfaprim * np.pi / 180)
print(f'Fx={Fx}, Fy={Fy}')

# Funkcja do obliczania reakcji podporowych
def reakcje(x1, x2, x3, d, Q, F0):
    x2prim = x2 - x1
    x3prim = x3 - x1
    dprim = d - x1
    R2 = (Q * (dprim - x2prim) - x2prim * F0) / (x3prim - x2prim)
    R1 = Q + F0 - R2
    return R1, R2

# Oblicz reakcje
R1x, R2x = reakcje(x1, x2, x3, d, 0, Fx)  # bo na kierunku x Q=0
R1y, R2y = reakcje(x1, x2, x3, d, Q, Fy)
print(f"R1x = {R1x:.2f} N, R2x = {R2x:.2f} N")
print(f"R1y = {R1y:.2f} N, R2y = {R2y:.2f} N")

# Generowanie punktów x
x = np.linspace(0, a, 500)

# Funkcja siły tnącej
def shear_force(x, R1, R2, Q, F0, x1, x2, d):
    T = np.zeros_like(x)
    for i, xi in enumerate(x):
        if xi < x1:
            T[i] = 0
        if xi >= x1:
            T[i] = F0
        if xi >= x2:
            T[i] = F0 - R1
        if xi >= x3:
            T[i] = F0 - R2 - R1
        if xi > d:
            T[i] = 0
    return T

# Funkcja momentu gnącego
def bending_moment(x, R1, R2, Q, F0, x1, x2, x3, d,T):
    M = np.zeros_like(x)
    for i, xi in enumerate(x):
        if xi < x1:
            M[i] = 0
        if xi >= x1:
            M[i] = F0*x1-xi*T[i]
        if xi >= x2:
            M[i] = F0*x1-T[i]*xi-R1*x2
        if xi >= x3:
            M[i] =  F0*x1-T[i]*xi-R1*x2-R2*x3
        if xi > d:
            M[i] = 0
    return M

def shearing_moment(x, R1, R2, Q, F0, x1, x2, x3, d,T):
    M=np.zeros_like(x)
    for i,xi in enumerate(x):
        if xi > x1 and xi < d:
            M[i] = M_sil*6/7 *0.94
        else:
            M[i]=0
    return M
# Oblicza siłę tnącą i moment gnący
Tx = shear_force(x, R1x, R2x, 0, Fx, x1, x2, d)
Mx = bending_moment(x, R1x, R2x, 0, Fx, x1, x2, x3, d, Tx)
M_s=shearing_moment(x, R1x, R2x, 0, F_0, x1, x2, x3, d, Tx)
Ty = shear_force(x, R1y, R2y, Q, Fy, x1, x2, d)
My = bending_moment(x, R1y, R2y, Q, Fy, x1, x2,x3, d, Ty)

# Sprawdzam moment na końcu wału (powinien być 0)
print(f"Moment na końcu wału (x=a): {np.sqrt(Mx[-1]**2 + My[-1]**2)} Nm")

T_wyp = [np.sqrt(Tx[i]**2 + Ty[i]**2) for i in range(len(x))]
M_wyp_g = [np.sqrt(Mx[i]**2 + My[i]**2) for i in range(len(x))]
M_wyp=[np.sqrt(M_wyp_g[i]**2 + (np.sqrt(3)*-.25*M_s[i])**2) for i in range(len(x))]

# Obliczanie wymaganej średnicy wału
tau_dopuszczalne = 0.1 * Zg * 1e6  # [Pa]
d_wymagany = [((32 * M_wyp[i]) / (np.pi * tau_dopuszczalne))**(1/3) * 1000 for i in range(len(x))]  # [mm]

# Wykresy (przystosowane do druku niskiej jakości)
mpl.rcParams.update({'font.size': 12, 'lines.linewidth': 1.2, 'axes.labelsize': 14, 'axes.titlesize': 14})

plt.figure(figsize=(10, 7))
plt.subplot(2, 1, 1)
plt.plot(x, Tx, color='black', linestyle='-', label='Tx(x) [N]')
plt.plot(x, Ty, color='grey', linestyle='--', label='Ty(x) [N]')
plt.plot(x, T_wyp, color='black', linestyle=':', label='T_wyp(x) [N]')
plt.axhline(0, color='black', linestyle=':')
plt.grid(True, color='grey', linestyle=':', linewidth=0.5)
plt.legend(loc='best', frameon=False)
plt.ylabel('Siła tnąca [N]')

plt.subplot(2, 1, 2)
plt.plot(x, Mx, color='black', linestyle='-', label='Mx(x) [Nm]')
plt.plot(x, My, color='grey', linestyle='--', label='My(x) [Nm]')
plt.plot(x, M_wyp, color='black', linestyle=':', label='M_wyp(x) [Nm]')
plt.plot(x, M_s, color='black', linestyle='-.', label='Ms(x) [Nm]')
plt.axhline(0, color='black', linestyle=':')
plt.legend(loc='best', frameon=False)
plt.grid(True, color='grey', linestyle=':', linewidth=0.5)
plt.ylabel('Moment [Nm]')
plt.xlabel('Położenie na osi wału [m]')
plt.tight_layout()
plt.savefig("M_T_x_print.png", dpi=200, bbox_inches='tight', facecolor='white')
plt.show()

plt.figure(figsize=(8, 4))
plt.plot(x, M_wyp, color='black', linestyle='-', label='M_wyp(x) [Nm]')
plt.xlabel('Położenie na osi wału [m]')
plt.ylabel('Wypadkowy moment gnący [Nm]')
plt.title('Rozkład wypadkowego momentu gnącego')
plt.grid(True, color='grey', linestyle=':', linewidth=0.5)
plt.legend(loc='best', frameon=False)
plt.tight_layout()
plt.savefig('Moment_wypadkowy_print.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.show()

# Definiowanie geometrii wału (x w metrach, średnice w mm)
geometria_walu = [
    (0.0, 28),
    (0.0625, 35),
    (0.095, 45),
    (0.33-0.105, 35),
    (0.33-0.06, 30),
    (0.330, 30)
]

# Generowanie obrysu wału
d_rzeczywisty = np.zeros_like(x)
current_d = geometria_walu[0][1]
for i, xi in enumerate(x):
    for point in geometria_walu:
        if xi >= point[0]:
            current_d = point[1]
    d_rzeczywisty[i] = current_d

plt.figure(figsize=(12, 6))
plt.plot(x, d_wymagany, 'b', label='Wymagana średnica (teoretyczna)')
plt.plot(x, d_rzeczywisty, 'r--', label='wartość średnic rzeczywistego wału', linewidth=2)
plt.ylabel('Średnica [mm]')
plt.xlabel('Położenie na osi wału [m]')
plt.axhline(0, color='k', linestyle='-')
plt.grid(True)
plt.legend()
plt.savefig('Srednica.png')
plt.show()

# Funkcja obliczająca naprężenia
def compute_stress(x, M_wyp, T_wyp, d_rzeczywisty, Zg, T_skret):
    """
    Oblicza naprężenia zginające, ścinające, skręcające i zastępcze (von Misesa).
    
    Parametry:
    T_skret : float
        Moment skręcający [Nm]
    """
    sigma = np.zeros_like(x)
    tau_transverse = np.zeros_like(x)
    tau_skret = np.zeros_like(x)
    von_mises = np.zeros_like(x)
    
    for i in range(len(x)):
        d = d_rzeczywisty[i] / 1000  # Konwersja mm -> m
        if d <= 0:
            continue
            
        # Obliczenia dla przekroju kołowego
        A = np.pi * (d**2) / 4      # Pole przekroju [m²]
        I = np.pi * (d**4) / 64     # Moment bezwładności [m^4]
        c = d / 2                   # Odległość od osi do skrajnego włókna [m]
        
        # Naprężenia zginające
        sigma_pa = (M_wyp[i] * c) / I  # [Pa]
        
        # Naprężenia ścinające poprzeczne
        tau_transverse_pa = (4/3) * (T_wyp[i] / A)  # [Pa]
        
        # Naprężenia skręcające
        tau_skret_pa = (16 * T_skret) / (np.pi * d**3)  # [Pa]
        
        # Naprężenia zastępcze (von Misesa)
        sigma_vm_pa = np.sqrt(sigma_pa**2 + 3 * (tau_transverse_pa**2 + tau_skret_pa**2))
        
        # Konwersja na MPa
        sigma[i] = sigma_pa / 1e6
        tau_transverse[i] = tau_transverse_pa / 1e6
        tau_skret[i] = tau_skret_pa / 1e6
        von_mises[i] = sigma_vm_pa / 1e6
    
    return {
        'sigma': sigma,
        'tau_transverse': tau_transverse,
        'tau_skret': tau_skret,
        'von_mises': von_mises,
        'Zg': Zg * np.ones_like(x)
    }

# Oblicz naprężenia z uwzględnieniem momentu skręcającego
stresses = compute_stress(x, M_wyp, T_wyp, d_rzeczywisty, Zg, T_skret=M_sil)# Oblicz naprężenia

plt.figure(figsize=(12, 6))
plt.plot(x, stresses['sigma'], 'b-', label='Naprężenia zginające σ [MPa]')
plt.plot(x, stresses['tau_transverse'], 'g-', label='Naprężenia ścinające poprzeczne τ [MPa]')
plt.plot(x, stresses['tau_skret'], 'm-', label='Naprężenia skręcające τ [MPa]')
plt.plot(x, stresses['von_mises'], 'r--', label='Naprężenia zastępcze (von Misesa) [MPa]')
plt.plot(x, stresses['Zg'], 'k:', label='Dopuszczalne naprężenia Zg [MPa]', linewidth=2)

plt.xlabel('Położenie na osi wału [m]')
plt.ylabel('Naprężenia [MPa]')
plt.title('Rozkład naprężeń w wale z uwzględnieniem skręcania')
plt.grid(True)
plt.legend()
plt.show()
print(f"Moment skręcający na wale: {M_sil*6/7*0.94} Nm")
R_1=np.sqrt(R1x**2+R1y**2)  # reakcja w pierwszej podporze
print(f"Reakcja w pierwszej podporze: {R_1:.2f} N")


#prędkość krytyczna wału
a=75/1000
b=162/1000
J=np.pi/64*((35/1000)**4)
E=210000 
f=(Q*a**2)/((3*E*J)*(a+3*b/2))/1000
print(f"Strzałka ugięcia: {f:.2f} mm")
omega_krytyczna = np.sqrt(9.81) / (f*10**-3)  # rad/s
print(f"Prędkość krytyczna wału: {omega_krytyczna:.2f} rad/s")
print(f"Prędkość krytyczna wału: {omega_krytyczna * 60 / (2 * np.pi):.2f} obr/min")