#include "DEInteg.h"
#include <cmath>
#include <stdexcept>
#include <limits>

/*
%----------------------------------------------------------------------------
%
% Purpose:
%   Numerical integration methods for ordinaray differential equations
%
%   This module provides implemenation of the variable order variable
%   stepsize multistep method of Shampine & Gordon.
%
% Last modified:   2015/08/25   M. Mahooti
%
% Reference:
%
%   Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
%   Freeman and Comp., San Francisco (1975).
%
%----------------------------------------------------------------------------*/

/*
Matrix DEInteg(
        const std::function<void(double, const Matrix &, Matrix &)> &func,
        double t0,
        const std::vector<double> &tout,
        double relerr,
        double abserr,
        const Matrix &y0
) {
    int n_eqn = y0.getFilas();
    int nt = static_cast<int>(tout.size());
    if (y0.getColumnas() != 1)
        throw std::invalid_argument("y0 debe ser n_eqn×1");
    // Preparar matriz de salida
    Matrix Yout(n_eqn, nt);

    // Parámetros de máquina
    const double twou = 2.0 * std::numeric_limits<double>::epsilon();
    const double fouru = 4.0 * std::numeric_limits<double>::epsilon();

    // Estados internos
    enum DE_STATE {
        DE_INIT = 1,
        DE_DONE = 2,
        DE_BADACC = 3,
        DE_NUMSTEPS = 4,
        DE_STIFF = 5,
        DE_INVPARAM = 6
    };

    int State = DE_INIT;
    bool PermitTOUT = true;
    double told = 0.0;

    // Potencias de dos
    std::vector<double> two = {
            1., 2., 4., 8., 16., 32., 64., 128.,
            256., 512., 1024., 2048., 4096., 8192.
    };
    // Constantes de la regla de error
    std::vector<double> gstr = {
            1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188,
            0.0143, 0.0114, 0.00936, 0.00789, 0.00679,
            0.00592, 0.00524, 0.00468
    };

    // Vectores de trabajo
    Matrix yy(n_eqn, 1), wt(n_eqn, 1), p(n_eqn, 1), yp(n_eqn, 1);
    Matrix phi(n_eqn, 17), g(14, 1), sig(14, 1), rho(14, 1),
            w(13, 1), alpha(13, 1), beta(13, 1), v(13, 1), psi(13, 1);

    // Copiamos y0 en yy
    for (int i = 1; i <= n_eqn; ++i) {
        yy(i, 1) = y0(i, 1);
    }

    // Índice de salida
    int k = 0;

    // Si t0 coincide con tout[0], devolvemos y0 directamente
    if (nt > 0 && std::abs(t0 - tout[0]) < 0.0) {
        Yout.setColumn(1, yy);
        k = 1;
    }

    // Inicializar paso h (por ejemplo h = (tout[nt-1]-t0)/10 o similar)
    double h = (nt > 0 ? (tout[0] - t0) : 0.0) / 10.0;
    if (h <= 0) h = 1.0;  // valor por defecto

    double t = t0;

    // --------------------------------------------------
    // Aquí iría el bucle principal while(k<nt && State==DE_INIT)
    //   - estimar yp = f(t, yy)
    //   - controlar errores con bob Check,
    //   - adaptar h,
    //   - avanzar t, yy,
    //   - almacenar en Yout(:,k)
    // --------------------------------------------------

    // Por ahora devolvemos ceros para compilar
    return Yout;
}

// Comprobación de parámetros inválidos
{
double epsilon = std::max(relerr, abserr);
if (
relerr<0.0 ||
       abserr < 0.0 ||
       epsilon <= 0.0 ||
       State> DE_INVPARAM
||
((State != DE_INIT) && (
std::abs(t
- told) > 0.0)))
{
State = DE_INVPARAM;
return
Yout;
}
}

// Definir intervalo de integración
double del = tout[k] - t;               // asumimos tout[k] es el próximo tiempo
double absdel = std::abs(del);
double tend = t + 100.0 * del;
if (!PermitTOUT) {
tend = tout[k];
}

// Inicializar contadores
int nostep = 0;
int kle4 = 0;
bool stiff = false;
double releps = relerr / epsilon;
double abseps = abserr / epsilon;

// Si arrancamos o la dirección cambió
if (State == DE_INIT || !OldPermit ||
delsgn *del
<= 0.0) {
start = true;
x = t;
yy = y0;                      // y0 es tu vector de estado actual
delsgn = (del >= 0.0 ? 1.0 : -1.0);
// paso inicial h = signo(max(4u*|x|,|tout-x|), tout-x)
double fouru = 4.0 * std::numeric_limits<double>::epsilon();
double h0 = std::max(fouru * std::abs(x), std::abs(tout[k] - x));
h = (tout[k] - x >= 0.0 ? h0 : -h0);
}

while (true) {
// Si hemos alcanzado o pasado el punto de salida
if (
std::abs(x
- tout[k]) >= absdel) {
// Inicializar acumuladores
Matrix yout(n_eqn, 1), ypout(n_eqn, 1);
g.assign(14, 0.0);
rho.assign(14, 0.0);
g[1]   = 1.0;  // MATLAB g(2)
rho[1] = 1.0;  // MATLAB rho(2)

double hi = tout[k] - x;
int ki = kold + 1;

// Inicializar w para el cómputo de g
for (
int i = 1;
i <=
ki;
++i) {
w[i] = 1.0 /
double(i);
}

// Calcular g[j] y rho[j]
double term = 0.0;
for (
int j = 1;
j <=
ki;
++j) {  // MATLAB j=2:ki
double psijm1 = psi_[j];
double gamma = (hi + term) / psijm1;
double eta = hi / psijm1;
for (
int i = 1;
i <= ki - j + 1; ++i) {
w[i] =
gamma *w[i]
-
eta *w[i + 1];
}
g[j+1]   = w[1];              // MATLAB g(j+1)=w(2)
rho[j+1] =
gamma *rho[j];    // rho(j+1)=gamma*rho(j)
term = psijm1;
}

// Interpolación de yout y ypout
for (
int j = 1;
j <=
ki;
++j) {
int i = ki + 1 - j;
// phi.col(i+1) está en phi.getColumn(i+1)
Matrix col = phi.getSubMatrix(1, n_eqn, i + 1, i + 1);
yout  +=
col *g[i + 1];
ypout +=
col *rho[i + 1];
}
// Ajustar por y + hi * sum(...)
yout = yy + yout * hi;
yy = yout;

// Estado final
State = DE_DONE;
x = tout[k];
told = x;
OldPermit = PermitTOUT;
return
yout_matrix; // o guarda yy en Yout(:,k)
}

// ... resto del bucle de pasos ...
}

// --- Extrapolación si no se permite pasar tout y estamos muy cerca ---
if (!
PermitTOUT &&std::abs(tout[k] - x)
<
fouru *std::abs(x)
) {
double hi = tout[k] - x;
Matrix yp(n_eqn, 1);
func(x, yy, yp
);               // derivada en x
y = yy + yp * hi;              // extrapola hasta tout
yy = y;
State = DE_DONE;
x = tout[k];
told = x;
OldPermit = PermitTOUT;
return
y;                      // o almacena en Yout y sigue
}

// --- (opcional) control de número máximo de pasos ---

// --- Limitar tamaño de paso y calcular vector de pesos ---
double remaining = tout[k] - x;
double hmax = std::min(std::abs(h), std::abs(remaining));
h = sign_(hmax, h);  // conserva signo de h
for (
int l = 1;
l <=
n_eqn;
++l) {
// si usas Matrix wt:
wt(l,
1) =
releps *std::abs(yy(l, 1))
+
abseps;
}

// --- Verificar paso demasiado pequeño para máquina ---
if (
std::abs(h)
<
fouru *std::abs(x)
) {
h = sign_(fouru * std::abs(x), h);
crash = true;
return

Matrix();  // salir con fallo
}

double p5eps = 0.5 * epsilon;
crash = false;
// Inicializar coeficientes g y sig para el primer step
g[1]   = 1.0;  // MATLAB g(2)
g[2]   = 0.5;  // MATLAB g(3)
sig[1] = 1.0;  // MATLAB sig(2)

int ifail = 0;

// ... aquí continúa el bloque "STEP" de tu integrador ...

// --- Evaluar si las tolerancias son demasiado estrictas ---
{
double round_ = 0.0;
for (
int l = 1;
l <=
n_eqn;
++l) {
double yy_l = yy(l, 1);
double w_l = wt(l - 1);  // si wt es vector 0-based
round_ += (
yy_l *yy_l
)/(
w_l *w_l
);
}
round_ = twou * std::sqrt(round_);
if (p5eps<round_) {
epsilon = 2.0 * round_ * (1.0 + fouru);
crash = true;
return
Yout;  // o lanzar excepción, según convención
}
}

// --- Bloque 0: inicialización del primer paso ---
if (start) {
// Evaluar derivada en el punto inicial
func(x, yy, yp
);

double sum = 0.0;
for (
int l = 1;
l <=
n_eqn;
++l) {
double yp_l = yp(l, 1);
double w_l = wt(l - 1);
phi(l,
2) =
yp_l;       // MATLAB phi(l,2)
phi(l,
3) = 0.0;        // phi(l,3)
sum += (
yp_l *yp_l
)/(
w_l *w_l
);
}
sum = std::sqrt(sum);

double absh = std::abs(h);
if (epsilon < 16.0 *
sum *h
* h) {
absh = 0.25 * std::sqrt(epsilon / sum);
}
h = sign_(std::max(absh, fouru * std::abs(x)), h);

// Inicializar contadores y flags
double hold = 0.0;
double hnew = 0.0;
int k = 1;
int kold = 0;
start = false;
phase1 = true;
bool nornd = true;
if (p5eps <= 100.0 * round_) {
nornd = false;
// setear phi(:,16)=0
for (
int l = 1;
l <=
n_eqn;
++l) {
phi(l,
16) = 0.0;
}
}
}
// --- Bloque 1: recalcular coeficientes si cambió h ---
{
int kp1 = k + 1;
int kp2 = k + 2;
int km1 = k - 1;
int km2 = k - 2;

// ns = número de pasos con tamaño h incluido el actual
if (h != hold) {
ns = 0;
}
if (ns <= kold) {
ns++;
}
int nsp1 = ns + 1;

if (k >= ns) {
// actualizar α, β, ψ, σ para los índices necesarios
beta[ns]     = 1.0;              // β(ns+1) en MATLAB
double realns = double(ns);
alpha[ns]    = 1.0/
realns;       // α(ns+1)
double temp1 = h * realns;
sig[nsp1]    = 1.0;              // σ(nsp1+1)

if (k >= nsp1) {
for (
int i = nsp1;
i <=
k;
++i) {
int im1 = i - 1;
double temp2 = psi_[im1];
psi_[im1]   =
temp1;                 // ψ(im1+1)
beta[i]     = beta[im1] * psi_[im1] /
temp2;  // β(i+1)
temp1 = temp2 + h;
alpha[i]    = h /
temp1;             // α(i+1)
double reali = double(i);
sig[i+1]    =
reali *alpha[i]
* sig[i]; // σ(i+2)
}
}
psi_[k] =
temp1;                         // ψ(k+1)

// calcular g[*]:
if (ns > 1) {
// si aumentó el orden, ajustar diagonal de v[*]
if (k > kold) {
double temp4 = double(k) * double(kp1);
v[k] = 1.0 /
temp4;               // v(k+1)
int nsm2 = ns - 2;
for (
int j = 1;
j <=
nsm2;
++j) {
int i = k - j;
v[i] = v[i] - alpha[j+1] * v[i+1];
}
}
// actualizar v[*] y w[*]
int limit1 = kp1 - ns;
double temp5 = alpha[ns];
for (
int iq = 1;
iq <=
limit1;
++iq) {
v[iq] = v[iq] -
temp5 *v[iq + 1];
w[iq] = v[iq];
}
g[nsp1] = w[1];  // g(nsp1+1)
}
else {
// orden inicial k ≤ 1
for (
int iq = 1;
iq <=
k;
++iq) {
double temp3 = double(iq * (iq + 1));
v[iq] = 1.0 /
temp3;
w[iq] = v[iq];
}
}

// terminar de calcular g[*]
int nsp2 = ns + 2;
if (kp1 >= nsp2) {
for (
int i = nsp2;
i <=
kp1;
++i) {
int limit2 = kp2 - i;
double temp6 = alpha[i - 1];
for (
int iq = 1;
iq <=
limit2;
++iq) {
w[iq] = w[iq] -
temp6 *w[iq + 1];
}
g[i] = w[1];  // g(i+1)
}
}
}
}
// fin Bloque 1
// === Begin Block 2 ===
// si k >= nsp1, actualizar phi star (multiplicar por beta)
if (k >= nsp1) {
for (
int i = nsp1;
i <=
k;
++i) {
double temp1 = beta[i];            // beta(i+1)
for (
int l = 1;
l <=
n_eqn;
++l) {
phi(l, i
+1) =
temp1 *phi(l, i + 1);
}
}
}

// Predecir solución p[] y limpiar phi para diferencias
for (
int l = 1;
l <=
n_eqn;
++l) {
phi(l, kp2
+1) =
phi(l, kp1
+1);
phi(l, kp1
+1) = 0.0;
p[l-1]        = 0.0;                // vector p es 0-based
}

// Combinar etapas para formar p
for (
int j = 1;
j <=
k;
++j) {
int i = kp1 - j;
int ip1 = i + 1;
double coeff = g[i];                // g(i+1)
for (
int l = 1;
l <=
n_eqn;
++l) {
p[l-1]     +=
coeff *phi(l, i + 1);
phi(l, i
+1) =
phi(l, i
+1) +
phi(l, ip1
+1);
}
}

// Aplicar corrección según nornd
if (nornd) {
for (
int l = 1;
l <=
n_eqn;
++l) {
p[l-1] =
yy(l,
1) +
h *p[l - 1];
}
} else {
for (
int l = 1;
l <=
n_eqn;
++l) {
double tau = h * p[l - 1] - phi(l, 16);
double ynew = yy(l, 1) + tau;
// almacenar corrección de redondeo si necesitas:
phi(l,
17) = (ynew -
yy(l,
1)) -
tau;
p[l-1]    =
ynew;
}
}

// avanzar el tiempo y evaluar derivadas nuevas
double xold = x;
x +=
h;
double absh = std::abs(h);
Matrix yp_mat(n_eqn, 1);
func(x, Matrix::fromVector(p), yp_mat
); // necesitas un constructor o método
// reemplaza con tu forma de convertir p[] a Matrix yp

// Estimar errores erkm2, erkm1, erk
double erkm2 = 0.0, erkm1 = 0.0, erk = 0.0;
for (
int l = 1;
l <=
n_eqn;
++l) {
double temp3 = 1.0 / wt[l - 1];
double temp4 = yp_mat(l, 1) - phi(l, 2); // phi(l,1+1)
if (km2 > 0) {
double sum = phi(l, km1 + 1) + temp4;  // phi(l,km1+1)
erkm2 += (
sum *temp3
)*(
sum *temp3
);
}
if (km2 >= 0) {
double sum = phi(l, k + 1) + temp4;    // phi(l,k+1)
erkm1 += (
sum *temp3
)*(
sum *temp3
);
}
erk   += (
temp4 *temp3
)*(
temp4 *temp3
);
}
if (km2 > 0) {
erkm2 = absh * sig[km1] * gstr[km2] * std::sqrt(erkm2);
}
if (km2 >= 0) {
erkm1 = absh * sig[k] * gstr[km1] * std::sqrt(erkm1);
}
double temp5 = absh * std::sqrt(erk);
double err = temp5 * (g[k] - g[kp1]);   // g(k+1)-g(kp1+1)
erk = temp5 * sig[kp1] * gstr[k]; // erk a orden k
int knew = k;

// decidir si bajar orden
if (km2 > 0) {
if (
std::max(erkm1, erkm2
) <= erk) {
knew = km1;
}
} else if (erkm1 <= 0.5*erk) {
knew = km1;
}

// === End Block 2 ===
// === Bloque 2: predicción y estimación de error ===

// 1) Predicción p[] y actualización de phi
for (
int l = 1;
l <=
n_eqn;
++l) {
// phi(l,kp2+1) = phi(l,kp1+1);
phi(l, kp2
+1) =
phi(l, kp1
+1);
// phi(l,kp1+1) = 0.0;
phi(l, kp1
+1) = 0.0;
// p(l) = 0.0;   -> p es vector 0-based
p[l-1] = 0.0;
}

for (
int j = 1;
j <=
k;
++j) {
int i = kp1 - j;
int ip1 = i + 1;
double coeff = g[i];  // g(i+1) en MATLAB
for (
int l = 1;
l <=
n_eqn;
++l) {
// p(l)     += coeff * phi(l,i+1);
p[l-1] +=
coeff *phi(l, i + 1);
// phi(l,i+1) += phi(l,ip1+1);
phi(l, i
+1) +=
phi(l, ip1
+1);
}
}

// 2) Corrección final según nornd
if (nornd) {
for (
int l = 1;
l <=
n_eqn;
++l) {
// p(l) = y(l) + h * p(l);
p[l-1] =
yy(l,
1) +
h *p[l - 1];
}
} else {
for (
int l = 1;
l <=
n_eqn;
++l) {
double tau = h * p[l - 1] - phi(l, 16);
double ynew = yy(l, 1) + tau;
// p(l) = ynew
p[l-1] =
ynew;
// phi(l,17) = (p(l)-y(l)) - tau
phi(l,
17) = (ynew -
yy(l,
1)) -
tau;
}
}

// 3) Avanzar x y evaluar derivada
double xold = x;
x +=
h;
double absh = std::abs(h);

// Construir Matrix para yp
Matrix yp_mat(n_eqn, 1);
{
Matrix p_mat(n_eqn, 1);
for (
int l = 1;
l <=
n_eqn;
++l)
p_mat(l,
1) = p[l-1];
func(x, p_mat, yp_mat
);
}

// 4) Estimación de errores erkm2, erkm1, erk
double erkm2 = 0.0, erkm1 = 0.0, erk = 0.0;
for (
int l = 1;
l <=
n_eqn;
++l) {
double tmp3 = 1.0 / wt[l - 1];
double tmp4 = yp_mat(l, 1) - phi(l, 2);  // phi(l,1+1)
if (km2 > 0) {
double sum = phi(l, km1 + 1) + tmp4;
erkm2 += (
sum *tmp3
)*(
sum *tmp3
);
}
if (km2 >= 0) {
double sum = phi(l, k + 1) + tmp4;
erkm1 += (
sum *tmp3
)*(
sum *tmp3
);
}
erk += (
tmp4 *tmp3
)*(
tmp4 *tmp3
);
}

if (km2 > 0) {
erkm2 = absh * sig[km1] * gstr[km2] * std::sqrt(erkm2);
}
if (km2 >= 0) {
erkm1 = absh * sig[k] * gstr[km1] * std::sqrt(erkm1);
}

double tmp5 = absh * std::sqrt(erk);
double err = tmp5 * (g[k] - g[kp1]);           // g(k+1)-g(kp1+1)
double erk_ord = tmp5 * sig[kp1] * gstr[k];     // erk a orden k
int knew = k;

// Decidir si bajamos de orden
if (km2 > 0) {
if (
std::max(erkm1, erkm2
) <= erk_ord) {
knew = km1;
}
} else {
if (erkm1 <= 0.5 * erk_ord) {
knew = km1;
}
}

// A partir de aquí vendría el test `if(err <= epsilon) …` o el bloque 3 en caso contrario.
// === Begin Block 3: manejo de paso fallido ===
if (!success) {
phase1 = false;
// Restaurar tiempo
x = xold;

// Restaurar phi: phi(:,i+1) = (phi(:,i+1) - phi(:,i+2)) / beta(i+1)
for (
int i = 1;
i <=
k;
++i) {
double invb = 1.0 / beta[i];     // beta(i+1)
int ip1 = i + 1;
for (
int l = 1;
l <=
n_eqn;
++l) {
double v1 = phi(l, i + 1) - phi(l, ip1 + 1);
phi(l, i
+1) =
invb *v1;
}
}

// Ajustar psi_: psi_(i) = psi_(i+1) - h
if (k >= 2) {
for (
int i = 2;
i <=
k;
++i) {
psi_[i-1] = psi_[i] -
h;      // ψ(i) = ψ(i+1) - h
}
}

// Contador de fallos consecutivos
ifail++;
double temp2 = 0.5;
if (ifail > 3) {
if (p5eps < 0.25 * erk) {
temp2 = std::sqrt(p5eps / erk);
}
}
if (ifail >= 3) {
knew = 1;   // bajar orden a 1
}

// Reducir paso
h = temp2 * h;
k = knew;
// Si el paso se vuelve demasiado pequeño
if (
std::abs(h)
<
fouru *std::abs(x)
) {
crash = true;
h = sign_(fouru * std::abs(x), h);
epsilon *= 2.0;
return
Yout;   // o lanza excepción / sale de la función
}
// Volver al inicio de Block 1
continue;
}
// === End Block 3 ===

// Si success, salimos del while para Block 4
if (success) {
break;
}
// === Begin Block 4 ===
// paso exitoso: actualizar orden y tamaño de paso para el siguiente
kold = k;
hold = h;

// corrección de la solución y evaluación de derivada
double tmp1 = h * g[kp1];  // g(kp1+1)
if (nornd) {
for (
int l = 1;
l <=
n_eqn;
++l) {
// y(l) = p(l) + tmp1*(yp(l) - phi(l,2));
y(l,
1) = p[l-1] +
tmp1 *(yp_mat(l, 1)
-
phi(l,
2));
}
} else {
for (
int l = 1;
l <=
n_eqn;
++l) {
double diff = yp_mat(l, 1) - phi(l, 2);
double rho = tmp1 * diff - phi(l, 17);
double ynew = p[l - 1] + rho;
y(l,
1)       =
ynew;
phi(l,
16)    = (ynew - p[l-1]) -
rho;
}
}
// reevaluar derivada en el punto corregido
func(x, y, yp_mat
);

// actualizar diferencias para el próximo paso
for (
int l = 1;
l <=
n_eqn;
++l) {
// phi(l,kp1+1) = yp(l) - phi(l,2);
phi(l, kp1
+1) =
yp_mat(l,
1) -
phi(l,
2);
// phi(l,kp2+1) = phi(l,kp1+1) - phi(l,kp2+1);
phi(l, kp2
+1) =
phi(l, kp1
+1) -
phi(l, kp2
+1);
}
// propagar corrección a todas las etapas
for (
int i = 1;
i <=
k;
++i) {
for (
int l = 1;
l <=
n_eqn;
++l) {
phi(l, i
+1) +=
phi(l, kp1
+1);
}
}

// estimar error de orden k+1 si procede
double erkp1 = 0.0;
if (knew == km1 || k == maxOrder) {
phase1 = false;
}
if (phase1) {
k = kp1;
erk = erkp1;
} else if (knew == km1) {
k = km1;
erk = erkm1;
} else if (kp1 <= ns) {
// calcular erkp1 = absh*gstr[kp1]*sqrt(sum((phi(:,kp2+1)/wt).^2))
for (
int l = 1;
l <=
n_eqn;
++l) {
double val = phi(l, kp2 + 1) / wt[l - 1];
erkp1 +=
val *val;
}
erkp1 = absh * gstr[kp1] * std::sqrt(erkp1);

if (k > 1) {
if (erkm1 <=
std::min(erk, erkp1
)) {
k = km1;
erk = erkm1;
} else if (erkp1<erk && k != maxOrder) {
k = kp1;
erk = erkp1;
}
} else if (erkp1 < 0.5 * erk) {
k = kp1;
erk = erkp1;
}
}
// === End Block 4 ===
// --- Actualizar h para el siguiente paso ---
if (phase1 || (p5eps >=
erk *two[k + 1]
)) {
hnew = 2.0 * h;
} else {
if (p5eps<erk) {
double temp2 = double(k + 1);
double r = std::pow(p5eps / erk, 1.0 / temp2);
double factor = std::max(0.5, std::min(0.9, r));
hnew = absh * factor;
hnew = sign_(std::max(hnew, fouru * std::abs(x)), h);
} else {
hnew = h;
}
}
h = hnew;

// --- Si hubo “crash” por tolerancias muy estrictas, salir con código de error ---
if (crash) {
State_ = DE_STATE::DE_BADACC;
relerr = epsilon * releps;
abserr = epsilon * abseps;
y = yy;    // restaurar último paso válido
t = x;
told = t;
OldPermit = true;
return
y;          // o lanza excepción, según convención de tu API
}

// contabilizar pasos y detectar rigidez
nostep++;
kle4++;
if (kold > 4)
kle4 = 0;
if (kle4 >= 50)
stiff = true;

// fin del while(true) — continúa al siguiente paso
*/