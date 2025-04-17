// dont do compute F, G
// convert helper functions folder too
// everything below qdot and all quaternion helper functions

#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589793
// WGS-84 Constants
#define EARTH_SEMI_MAJOR_AXIS 6378137.0    // a (meters)
#define EARTH_ECCENTRICITY 0.0818191908425 // e

double magI[3];

double sind(double phi)
{
    return sin(phi * M_PI / 180);
}

double cosd(double phi)
{
    return cos(phi * M_PI / 180);
}

double secd(double phi)
{
    return 1 / cosd(phi);
}

double sec(double phi)
{
    return 1 / cos(phi);
}

double tand(double phi)
{
    return sind(phi) / cosd(phi);
}

void compute_radii(double *R_phi, double *R_lamb, double *dR_phi_dphi, double *dR_lamb_dphi, double x[22])
{
    long a = 6378137;
    double b = 6356752.31425;
    double e2 = 1 - pow(b / a, 2);
    double phi = x[4];

    double num1 = a * (1 - e2);
    double den1 = pow((1 - e2 * pow(sin(phi * M_PI / 180), 2)), 3 / 2);
    double den2 = pow((1 - e2 * pow(sin(phi * M_PI / 180), 2)), 1 / 2);

    *R_phi = num1 / den1;
    *R_lamb = a / den2;

    double num3 = 3 * a * (1 - e2) * e2 * sin(phi * M_PI / 180) * cos(phi * M_PI);
    double den3 = pow(1 - e2 * pow(sin(phi * M_PI / 180), 2), 5 / 2);

    double num4 = a * e2 * sin(phi * M_PI / 180) * cos(phi * M_PI / 180);
    double den4 = pow(1 - e2 * pow(sin(phi * M_PI / 180), 2), 3 / 2);

    *dR_lamb_dphi = num4 / den4;
    *dR_phi_dphi = num3 / num3;
}

void compute_dwdp(double result[3][3], double x[22], double we)
{
    double phi = x[4];
    double h = x[6];
    double vn = x[7];
    double ve = x[8];

    double R_phi;
    double R_lamb;
    double dR_phi_dphi;
    double dR_lamb_dphi;

    compute_radii(&R_phi, &R_lamb, &dR_phi_dphi, &dR_lamb_dphi, x);

    double m11 = -we * sind(phi) - ve / pow(R_lamb + h, 2) * dR_lamb_dphi;
    double m13 = -ve / pow(R_lamb + h, 2);
    double m21 = vn / pow(R_phi + h, 2) * dR_phi_dphi;
    double m23 = vn / pow(R_phi + h, 2);
    double m31 = -we * cosd(phi) - (ve * pow(secd(phi), 2)) / (R_lamb + h) + ve * tand(phi) / pow(R_lamb + h, 2) * dR_lamb_dphi;
    double m33 = (ve * tand(phi)) / pow(R_lamb + h, 2);

    double dwdp[3][3] = {
        {m11, 0, m13},
        {m21, 0, m23},
        {m31, 0, m33}};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = dwdp[i][j];
        }
    }
}

void compute_dwdv(double result[3][3], double x[])
{
    double phi = x[4];
    double h = x[6];

    double R_phi;
    double R_lamb;
    double _a;
    double _b;

    compute_radii(&R_phi, &R_lamb, &_a, &_b, x);
    double m12 = 1 / (R_lamb + h);
    double m21 = -1 / (R_phi + h);
    double m32 = -tand(phi) / (R_lamb + h);

    double dwdv[3][3] = {
        {0, m12, 0},
        {m21, 0, 0},
        {0, m32, 0}};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = dwdv[i][j];
        }
    }
}

void compute_dpdot_dp(double result[3][3], double x[])
{
    double phi = x[4];
    double h = x[6];
    double vn = x[7];
    double ve = x[8];

    double R_phi;
    double R_lamb;
    double dR_phi_dphi;
    double dR_lamb_dphi;

    compute_radii(&R_phi, &R_lamb, &dR_phi_dphi, &dR_lamb_dphi, x);

    double m11 = -vn / pow(R_phi + h, 2) * dR_phi_dphi;
    double m13 = -vn / pow(R_phi + h, 2);
    double m21 = -(ve * sec(phi)) / pow(R_lamb + h, 2) * dR_lamb_dphi + ve * secd(phi) * tand(phi) / (R_lamb + h);
    double m23 = -ve * secd(phi) / pow(R_lamb + h, 2);

    double dpdot_dp[3][3] = {
        {m11, 0, m13},
        {m21, 0, m23},
        {0, 0, 0}};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = dpdot_dp[i][j];
        }
    }
}

void compute_dpdot_dv(double result[3][3], double x[22])
{
    double phi = x[4];
    double h = x[6];

    double R_phi;
    double R_lamb;
    double _a;
    double _b;

    compute_radii(&R_phi, &R_lamb, &_a, &_b, x);

    double m11 = 1 / (R_phi + h);
    double m22 = secd(phi) / (R_lamb + h);
    double m33 = -1;

    double dpdot_dv[3][3] = {
        {m11, 0.0, 0.0},
        {0.0, m22, 0.0},
        {0.0, 0.0, m33}};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = dpdot_dv[i][j];
        }
    }
}

void compute_g(double *g, double *dg_dphi, double *dg_dh, double x[22])
{
    double phi = x[4];
    double h = x[6];
    double grav = 9.780327 * (1 + 5.3024 * pow(10, -3) * pow(sind(phi), 2) - 5.8 * pow(10, -6) * pow(sind(2 * phi), 2)) - (3.0877 * pow(10, -6) - 4.4 * pow(10, -9) * pow(sind(phi), 2)) * h + 7.2 * pow(10, -14) * pow(h, 2);

    *dg_dphi = 9.780327 * (1.06048 * pow(10, -2) * sind(phi) * cosd(phi) - 4.64 * pow(10, -5) * (sind(phi) * pow(cosd(phi), 3) - pow(sind(phi), 3) * cosd(phi))) + 8.8 * pow(10, -9) * h * sind(phi) * cosd(phi);

    *dg_dh = -3.0877 * pow(10, -6) + 4.4 * pow(10, -9) * pow(sind(phi), 2) + 1.44 * pow(10, -13) * h;
}

void compute_dvdot_dp(double result[3][3], double x[22], double we)
{
    double phi = x[4];
    double h = x[6];
    double vn = x[7];
    double ve = x[8];
    double vd = x[9];

    double R_phi;
    double R_lamb;
    double dR_phi_dphi;
    double dR_lamb_dphi;

    compute_radii(&R_phi, &R_lamb, &dR_phi_dphi, &dR_lamb_dphi, x);

    double _c;
    double dg_dphi;
    double dg_dh;

    compute_g(&_c, &dg_dphi, &dg_dh, x);

    double Y11 = -(pow(ve, 2) * pow(sind(phi), 2)) / (R_lamb + h) + (pow(ve, 2) * tand(phi)) / pow(R_lamb + h, 2) * dR_lamb_dphi - 2 * we * ve * cosd(phi) - vn * vd / pow(R_phi + h, 2) * dR_phi_dphi;
    double Y13 = pow(ve, 2) * tand(phi) / pow(R_lamb + h, 2) - vn * vd / pow(R_phi + h, 2);
    double Y21 = ve * vn * pow(secd(phi), 2) / (R_lamb + h) - ve * vn * tand(phi) / pow(R_lamb + h, 2) * dR_lamb_dphi + 2 * we * vn * cosd(phi) - ve * vd / pow(R_lamb + h, 2) * dR_lamb_dphi - 2 * we * vd * sind(phi);
    double Y23 = -ve * ((vn * tand(phi) + vd) / pow(R_lamb + h, 2));
    double Y31 = pow(ve, 2) / pow(R_lamb + h, 2) * dR_lamb_dphi + pow(vn, 2) / pow(R_phi + h, 2) * dR_phi_dphi + 2 * we * ve * sind(phi) + dg_dphi;
    double Y33 = pow(ve, 2) / pow(R_lamb + h, 2) + pow(vn, 2) / pow(R_phi + h, 2) + dg_dh;

    double dvdot_dp[3][3] = {
        {Y11, 0, Y13},
        {Y21, 0, Y23},
        {Y31, 0, Y33}};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = dvdot_dp[i][j];
        }
    }
}

void compute_dvdot_dv(double result[3][3], double x[22], double we)
{
    double phi = x[4];
    double h = x[6];
    double vn = x[7];
    double ve = x[8];
    double vd = x[9];

    double R_phi;
    double R_lamb;
    double dR_phi_dphi;
    double dR_lamb_dphi;

    compute_radii(&R_phi, &R_lamb, &dR_phi_dphi, &dR_lamb_dphi, x);

    double Z11 = vd / (R_phi + h);
    double Z12 = (-2 * ve * tand(phi)) / (R_lamb + h) + 2 * we * sind(phi);
    double Z13 = vn / (R_phi + h);
    double Z21 = (ve * tand(phi)) / (R_lamb + h) + 2 * we * sind(phi);
    double Z22 = (vd + vn * tand(phi)) / (R_lamb + h);
    double Z23 = ve / (R_lamb + h) + 2 * we * cosd(phi);
    double Z31 = (-2 * vn) / (R_phi + h);
    double Z32 = (-2 * ve) / (R_lamb + h) - 2 * we * cosd(phi);

    double dvdot_dv[3][3] = {
        {Z11, Z12, Z13},
        {Z21, Z22, Z23},
        {Z31, Z32, 0}};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = dvdot_dv[i][j];
        }
    }
}

void createSkewSymmMat(double result[3][3], double vec[3])
{
    double ssMat[3][3] = {
        {0, -vec[2], vec[1]},
        {vec[2], 0, -vec[1]},
        {-vec[1], vec[0], 0}};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = result[i][j];
        }
    }
}

void integrate(double x[22], double *P, double *Pq, double dt,
               double *qdot, double *pdot, double *vdot,
               double *Pdot, double *Pqdot,
               double *x_minus, double *P_minus, double *Pq_minus)
{

    // Assemble xdot vector
    double xdot[16];
    for (int i = 0; i < 4; i++)
    {
        xdot[i] = qdot[i];
        xdot[i + 4] = pdot[i];
        xdot[i + 8] = vdot[i];
    }
    for (int i = 12; i < 16; i++)
    {
        xdot[i] = 0.0; // Initialize last 4 elements to zero
    }

    // Euler integration
    for (int i = 0; i < 16; i++)
    {
        x_minus[i] = x[i] + dt * xdot[i];
    }

    // Normalize quaternion (first 4 elements of x_minus)
    double q_norm = 0.0;
    for (int i = 0; i < 4; i++)
    {
        q_norm += x_minus[i] * x_minus[i];
    }
    q_norm = sqrt(q_norm);

    if (q_norm > 1.02)
    {
        // Warning message in case of large quaternion norm
        printf("Warning: Quaternion norm of %.2f is > 2 percent of unity, consider decreasing integrator step size.\n", q_norm);
    }

    // Normalize quaternion
    for (int i = 0; i < 4; i++)
    {
        x_minus[i] = x_minus[i] / q_norm;
    }

    // Euler integration for P and Pq
    for (int i = 0; i < 9; i++)
    {
        P_minus[i] = P[i] + dt * Pdot[i];
    }

    for (int i = 0; i < 9; i++)
    {
        Pq_minus[i] = Pq[i] + dt * Pqdot[i];
    }
}

void createXiMat(double result[4][3], double x[22])
{
    double q1 = x[0];
    double q2 = x[1];
    double q3 = x[2];
    double q4 = x[3];

    double Xi[4][3] = {
        {q1, -q4, q3},
        {q4, q1, -q2},
        {-q3, q2, q1},
        {-q2, -q3, -q4}};

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = Xi[i][j];
        }
    }
}

void compute_wn(double result[3], double x[22], double we)
{
    double phi = x[4];
    double h = x[6];
    double vn = x[7];
    double ve = x[8];

    double R_phi, R_lamb;
    double _x;
    double _y;
    compute_radii(&R_phi, &R_lamb, &_x, &_y, x);

    double vec1[3] = {we * cos(phi * M_PI / 180), 0, -we * sin(phi * M_PI / 180)};
    double vec2[3] = {ve / (R_lamb + h),
                      -vn / (R_phi + h),
                      -(ve * tan(phi * M_PI / 180) / (R_lamb + h))};

    double wn[3] = {vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]};
    result[0] = wn[0];
    result[1] = wn[1];
    result[2] = wn[2];
}

void quaternion2DCM(double DCM[3][3], double q[4])
{
    double qw = q[0], qx = q[1], qy = q[2], qz = q[3];

    // Compute the Direction Cosine Matrix (DCM)
    DCM[0][0] = 1 - 2 * (qy * qy + qz * qz);
    DCM[0][1] = 2 * (qx * qy - qz * qw);
    DCM[0][2] = 2 * (qx * qz + qy * qw);

    DCM[1][0] = 2 * (qx * qy + qz * qw);
    DCM[1][1] = 1 - 2 * (qx * qx + qz * qz);
    DCM[1][2] = 2 * (qy * qz - qx * qw);

    DCM[2][0] = 2 * (qx * qz - qy * qw);
    DCM[2][1] = 2 * (qy * qz + qx * qw);
    DCM[2][2] = 1 - 2 * (qx * qx + qy * qy);
}

void transpose_matrix(double result[3][3], double matrix[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[j][i] = matrix[i][j];
        }
    }
}

void compute_what(double result[3], double x[22], double w_meas[4], double we)
{
    double sf_g[3] = {0};
    for (int i = 0; i < 3; i++)
    {
        sf_g[i] = x[i + 16];
    }

    double bias_g[3] = {0};
    for (int i = 0; i < 3; i++)
    {
        bias_g[i] = x[10 + i];
    }

    double q[4] = {0};
    for (int i = 0; i < 4; i++)
    {
        q[i] = x[i];
    }

    double DCM[3][3];
    double D_bn[3][3];

    quaternion2DCM(DCM, q);
    transpose_matrix(D_bn, DCM);

    double w_n[3];
    compute_wn(w_n, x, we);

    double prod[3];
    for (int i = 0; i < 3; i++)
    {
        prod[i] = D_bn[i][0] * w_n[0] + D_bn[i][1] * w_n[1] + D_bn[i][2] * w_n[2];
    }

    double w_hat[3];
    for (int i = 0; i < 3; i++)
    {
        w_hat[i] = 1 / (1 + sf_g[i]) * (w_meas[i] - bias_g[i]) - prod[i];
    }

    for (int i = 0; i < 3; i++)
    {
        result[i] = w_hat[i];
    }
}

void copy_matrix(int size, double dest[size][size], double src[size][size])
{
    memcpy(dest, src, size * size * sizeof(double));
}

void copy_array(double *dest, const double *src, int size)
{
    memcpy(dest, src, size * sizeof(double));
}

void skew(double v[3], double S[3][3])
{
    S[0][0] = 0;
    S[0][1] = -v[2];
    S[0][2] = v[1];
    S[1][0] = v[2];
    S[1][1] = 0;
    S[1][2] = -v[0];
    S[2][0] = -v[1];
    S[2][1] = v[0];
    S[2][2] = 0;
}

double deg2rad(double deg)
{
    return deg * PI / 180;
}

// QUATERNIONS
double dot(double v1[3], double v2[3])
{
    double dot_product = 0;
    for (int i = 0; i < 3; i++)
    {
        dot_product += v1[i] * v2[i];
    }

    return dot_product;
}

void cross(double v1[3], double v2[3], double v_cross[3])
{
    v_cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v_cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v_cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void q2v(double q[4], double v[3])
{
    if (q[0] > 1e-6)
    {
        printf("q2v should only be used on pure quaternions.");
    }
    v[0] = q[1];
    v[1] = q[2];
    v[2] = q[3];
}

double norm(double q[4])
{
    double norm = 0.0;
    for (int i = 0; i < 4; i++)
    {
        norm += q[i] * q[i];
    }

    return norm;
}

void qconj(double q[4], double qbar[4])
{
    double norm_val = norm(q);
    norm_val = sqrt(norm_val);
    if (fabs(norm_val - 1) > 1e-6)
    {
        printf("Warning: Quaternion is not unit norm. Conjugate is not inverse.\n");
    }

    qbar[0] = q[0];
    for (int i = 1; i < 4; i++)
    {
        qbar[i] = -q[i];
    }
}

void qmult(double q1[4], double q2[4], double q_out[4])
{
    double s1 = q1[0];
    double v1[3];
    q2v(q1, v1);

    double s2 = q2[0];
    double v2[3];
    q2v(q2, v2);

    double s_out = s1 * s2 - dot(v1, v2);
    double v_cross[3];
    cross(q1, q2, v_cross);

    double prod1[3];
    for (int i = 0; i < 3; i++)
    {
        prod1[i] = s1 * v2[i];
    }
    double prod2[3];
    for (int i = 0; i < 3; i++)
    {
        prod1[i] = s2 * v1[i];
    }

    double v_out[3];

    for (int i = 0; i < 3; i++)
    {
        v_out[i] = prod1[i] + prod2[i] + v_cross[i];
    }

    q_out[0] = s_out;
    for (int i = 0; i < 3; i++)
    {
        q_out[i + 1] = v_out[i];
    }
}

void quaternion_exp(double v[3], double dq[4])
{
    double vnorm = 0;
    for (int i = 0; i < 3; i++)
    {
        vnorm += v[i] * v[i];
    }

    dq[0] = cos(vnorm);
    double s = sin(vnorm) / vnorm;
    for (int i = 0; i < 3; i++)
    {
        dq[i + 1] = s * v[i];
    }
}

void quaternion2DCM_new(double q[4], double CB21[3][3])
{
    double norm_q = norm(q);
    norm_q = sqrt(norm_q);

    for (int i = 0; i < 4; i++)
    {
        q[i] /= norm_q;
    }

    double s = q[0];
    double v[3] = {q[1], q[2], q[3]};

    double v_dot_v = 0.0;
    for (int i = 0; i < 3; i++)
    {
        v_dot_v += v[i] * v[i];
    }

    double S[3][3];
    skew(v, S);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            CB21[i][j] = (s * s - v_dot_v) * (i == j) + 2 * v[i] * v[j] + 2 * S[i][j];
        }
    }
}

void sandwich(double q[4], double x[22], double y[4])
{
    double cur[4];
    double qbar[4];
    qconj(q, qbar);
    qmult(x, qbar, cur);
    qmult(q, cur, y);
}

// HELPER FUNCTIONS
void ecef2lla(double ecef[3], double lla[3])
{
    double a = 6378137;
    double b = 6356752.3124;
    double e = 0.08181919;
    double e2 = e * e;
    double f = 1 / 298.257223563;
    double x = ecef[0];
    double y = ecef[1];
    double z = ecef[2];
    double p = x * x + y * y;
    double psi = atan2(z * a, p * b);
    double ebar2 = pow(a / b, 2) - 1;
    double phi = 180 / PI * atan2(z + ebar2 * b * pow(sin(psi), 3), p - e2 * a * pow(cos(psi), 3));

    double N = a / sqrt(1 - e2 * pow(sind(phi), 2));
    double lambda = 180 / PI * atan2(y, x);
    double h = p / cosd(phi) - N;

    lla[0] = phi;
    lla[1] = lambda;
    lla[2] = h;
}

void lla2ecef(double lla[3], double ecef[3])
{
    double a = 6378137;
    double b = 6356752.3124;
    double e = 0.08181919;
    double phi = lla[0];
    double lambda = lla[1];
    double h = lla[2];
    double N = a / sqrt(1 - pow(e * sind(phi), 2));
    double x = (N + h) * cosd(phi) * cosd(lambda);
    double y = (N + h) * cosd(phi) * sind(lambda);
    double z = (N * (1 - e * e) + h) * sind(phi);
    ecef[0] = x;
    ecef[1] = y;
    ecef[2] = z;
}

void ecef2ned(double ecef[3], double lla0[3], double ned[3])
{
    double ecef0[3];
    lla2ecef(lla0, ecef0);
    double phi = lla0[0];
    double lambda = lla0[1];
    double s1 = sind(phi);
    double c1 = cosd(phi);
    double s2 = sind(lambda);
    double c2 = cosd(lambda);
    double Ae2n[3][3] = {
        {-s1 * c2, -s1 * s2, c1}, {-s2, c2, 0}, {-c1 * c2, -c1 * s2, -s1}};

    for (int i = 0; i < 3; i++)
    {
        ned[i] = 0;
        for (int j = 0; j < 3; j++)
        {
            ned[i] += Ae2n[i][j] * (ecef[j] - ecef0[j]);
        }
    }
}

void improvedecef2lla(double r[3], double geo_lla[3])
{
    double a = 6378136.3;
    double f = 1 / 298.257;
    double b = a * (1 - f);
    double e = sqrt(1.0 - (b * b / a * a));
    double epsilon = pow(10, -12);
    int max_iter = 100;
    double x = r[0];
    double y = r[1];
    double z = r[2];
    double p = sqrt(x * x + y * y);
    double z_norm = fabs(z);
    double phi;
    double lambda;
    double h_geo;
    if (p == 0)
    {
        if (z == 0)
        {
            phi = 0.0;
            lambda = 0.0;
            h_geo = -b;
        }
        else
        {
            if (z > 0)
            {
                phi = PI / 2.0;
            }
            else
            {
                phi = -PI / 2.0;
            }
            lambda = 0.0;
            h_geo = z_norm - b;
        }
    }
    else if (z == 0)
    {
        phi = 0.0;
        lambda = atan2(y, x);
        h_geo = p - a;
    }
    else
    {
        double c = a * e * e;
        double e_prime = sqrt(1 - e * e);
        double z_prime = e_prime * z_norm;
        double u = 2.0 * (z_prime - c);
        double v = 2.0 * (z_prime + c);
        double t_M = (c - z_prime) / p;
        double f_M = p * pow(t_M, 4) + u * pow(t_M, 3) + v * t_M - p;
        double t_1 = (p - c + z_prime) / (p - c + 2.0 * z_prime);
        double t_0 = p / (z_prime + c);
        double t;

        if (t_M <= 0)
            t = t_1;
        else if (t_M >= 1)
            t = t_0;
        else if (f_M >= 0)
            t = t_0;
        else
            t = t_1;

        for (int i = 0; i < max_iter; i++)
        {
            double delta_t = (p - (p * pow(t, 4) + u * pow(t, 3) + v * t)) / (4.0 * p * pow(t, 3) + 3.0 * u * pow(t, 2) + v);
            double prev_t = t;
            t = t + delta_t;
            double TE = fabs(t - prev_t);
            if (TE <= epsilon)
                break;
        }

        phi = atan2(1.0 - pow(t, 2), 2.0 * e_prime * t);
        if (z <= 0)
        {
            phi = -phi;
        }

        h_geo = ((2.0 * p * e_prime * t) + (z_norm * (1.0 - t * t)) - (a * e_prime * (1.0 + t * t))) / (sqrt(pow(1.0 + t * t, 2) - (4.0 * e * e * t * t)));
        lambda = atan2(y, x);
        if (lambda < 0)
            lambda = 2.0 * PI + lambda;
    }
    geo_lla[0] = phi;
    geo_lla[1] = lambda;
    geo_lla[2] = h_geo;
}

void lla2ned(double lla[3], double lla0[3], double ned[3])
{
    double ecef[3];
    lla2ecef(lla, ecef);
    double result[3];
    ecef2ned(ecef, lla0, result);
    for (int i = 0; i < 3; i++)
    {
        ned[i] = result[i];
    }
}

void ned2ecef(double ned[3], double lla0[3], double ecef[3])
{
    double ecef0[3];
    lla2ecef(lla0, ecef0);
    double phi = lla0[0];
    double lambda = lla0[1];
    double s1 = sind(phi);
    double c1 = cosd(phi);
    double s2 = sind(lambda);
    double c2 = cosd(lambda);
    double An2e[3][3] = {
        {-s1 * c2, -s2, -c1 * c2},
        {-s1 * s2, c2, -c1 * s2},
        {c1, 0, -s1}};

    double en[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++)
    {
        en[i] = 0;
        for (int j = 0; j < 3; j++)
        {
            en[i] += An2e[i][j] * ned[j];
        }
    }

    for (int i = 0; i < 3; i++)
    {
        ecef[i] = en[i] + ecef0[i];
    }
}

void ned2lla(double ned[3], double lla0[3], double lla[3])
{
    double ecef[3];
    ned2ecef(ned, lla0, ecef);
    ecef2lla(ecef, lla);
}

int main()
{
    // DEFINE CONSANTS
    double we = 7.29211 * pow(10, -5); // Earth's sidereal rotation rate [rad/s]

    // DEFINE PARAMETERS
    int IMU_freq = 100;
    double dt = 1 / IMU_freq;
    int GPS_freq = 20;
    int mag_freq = 100;

    // define H
    double H[3][21] = {0};
    for (int i = 0; i < 3; i++)
    {
        H[i][3 + i] = 1.0;
    }

    double Hb = 0; // barometer measurement function Jacobian
    double Hq[3][6] = {0};
    double S[3][3];
    skew(magI, S);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Hq[i][j] = S[i][j];
        }
    }

    double R[3][3] = {
        {1.35e-5, 0.0, 0.0},
        {0.0, 1.65e-5, 0.0},
        {0.0, 0.0, 2.0}};

    double Rq[3][3] = {
        {3.2e-7, 0.0, 0.0},
        {0.0, 4.1e-7, 0.0},
        {0.0, 0.0, 3.2e-7}};

    double Rb = 1e-8;

    double scale = deg2rad(12e-3) * sqrt(IMU_freq);
    double nu_gv_mat[3][3] = {
        {scale, 0.0, 0.0},
        {0.0, scale, 0.0},
        {0.0, 0.0, scale}};

    scale = 1 / IMU_freq * deg2rad(3 / 3600);
    double nu_gu_mat[3][3] = {
        {scale, 0.0, 0.0},
        {0.0, scale, 0.0},
        {0.0, 0.0, scale}};

    scale = (200e-6 * 9.81) * sqrt(IMU_freq);
    double nu_av_mat[3][3] = {
        {scale, 0.0, 0.0},
        {0.0, scale, 0.0},
        {0.0, 0.0, scale}};

    scale = 1 / IMU_freq * 40e-6 * 9.8;
    double nu_au_mat[3][3] = {
        {scale, 0.0, 0.0},
        {0.0, scale, 0.0},
        {0.0, 0.0, scale}};

    double Q[12][12] = {0};
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Q[i][j] = nu_gv_mat[i][j];
        }
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Q[i + 3][j + 3] = nu_gu_mat[i][j];
        }
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Q[i + 6][j + 6] = nu_av_mat[i][j];
        }
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Q[i + 9][j + 9] = nu_au_mat[i][j];
        }
    }

    double Qq[6][6] = {0};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Qq[i][j] = nu_gv_mat[i][j];
        }
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Qq[i + 3][j + 3] = nu_gu_mat[i][j];
        }
    }

    double lla0[] = {30.9275, -81.51472222222, 45};
    double q0[] = {0.707106781186548, 0, 0.707106781186547, 0};

    double x0[19] = {0};
    for (int i = 0; i < 4; i++)
    {
        x0[i] = q0[i];
    }

    for (int i = 0; i < 3; i++)
    {
        x0[i + 4] = lla0[i];
    }

    double att_unc0 = 1e-4;   // initial attitude uncertainty [rad]
    double pos_unc0 = 1e-4;   // initial position uncertainty [m]
    double vel_unc0 = 1e-4;   // initial velocity uncertainty [m/s]
    double gbias_unc0 = 1e-4; // initial gyro bias uncertainty [rad/s]
    double abias_unc0 = 1e-4; // initial accel bias uncertainty [rad/s]
    double gsf_unc0 = 1e-4;   // gyro scale factor uncertainty [rad/s]
    double asf_unc0 = 1e-4;   // accel scale factor uncertainty [rad/s]

    double diag_values[21] = {
        att_unc0, att_unc0, att_unc0,
        pos_unc0, pos_unc0, pos_unc0,
        vel_unc0, vel_unc0, vel_unc0,
        gbias_unc0, gbias_unc0, gbias_unc0,
        abias_unc0, abias_unc0, abias_unc0,
        gsf_unc0, gsf_unc0, gsf_unc0,
        asf_unc0, asf_unc0, asf_unc0};

    double P0[19][19] = {0};
    for (int i = 0; i < 19; i++)
    {
        P0[i][i] = diag_values[i];
    }

    double pq0_diag_vals[6] = {
        att_unc0,
        att_unc0,
        att_unc0,
        abias_unc0,
        abias_unc0,
        abias_unc0};

    double Pq0[6][6] = {0};

    for (int i = 0; i < 6; i++)
    {
        Pq0[i][i] = diag_values[i];
    }

    double x_plus[19];
    copy_array(x_plus, x0, 19);

    double P_plus[18][18] = {0};
    copy_matrix(18, P_plus, P0);

    double Pq_plus[6][6] = {0};
    copy_matrix(6, Pq_plus, Pq0);

    int n, N;

    for (int i = 0; i < N; i++)
    {
        double a_meas;
        double w_meas;
    }

    return 0;
}
