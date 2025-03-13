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

void compute_radii(double *Re, double *Rn, double x[3])
{
    // Extract latitude (assumes x[1] contains latitude in radians)
    double lat = x[1];

    // Compute the prime vertical radius of curvature (Re)
    double sin_lat = sin(lat);
    *Re = EARTH_SEMI_MAJOR_AXIS / sqrt(1 - pow(EARTH_ECCENTRICITY, 2) * pow(sin_lat, 2));

    // Compute the meridional radius of curvature (Rn)
    *Rn = (*Re * (1 - pow(EARTH_ECCENTRICITY, 2))) / (1 - pow(EARTH_ECCENTRICITY, 2) * pow(sin_lat, 2));
}

void compute_wn(double x[], double we)
{
    double phi = x[4];
    double h = x[6];
    double vn = x[7];
    double ve = x[8];

    double R_phi, R_lamb;
    compute_radii(&R_phi, &R_lamb, x);

    double vec1[3] = {we * cos(phi * M_PI / 180), 0, -we * sin(phi * M_PI / 180)};
    double vec2[3] = {ve / (R_lamb + h),
                      -vn / (R_phi + h),
                      -(ve * tan(phi * M_PI / 180) / (R_lamb + h))};

    double wn[3] = vec1[0] + vec2[0] + vec1[1] + vec2[1] + vec1[2] + vec2[2];
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

double compute_what(double x[], double w_meas[], double we)
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

    double diag_values[19] = {
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
