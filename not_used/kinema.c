#include <stdio.h>
#include <math.h>

#define sqrt3 1.7320508075688772935274463415059
#define pi 3.1415926535897932384626433832795
#define sin120 sqrt3 / 2.0
#define cos120 -0.5
#define tan60 sqrt3
#define sin30 0.5
#define tan30 1.0 / sqrt3

#define non_existing_povar_error -2
#define no_error 1

typedef struct {
    double x, y, z;
    double a, b, c;
    double ArmLength, RodLength, BassTri, PlatformTri;
} DeltaKinematics;

void DeltaKinematics_setup(DeltaKinematics *dk, double _ArmLength, double _RodLength, double _BassTri, double _PlatformTri);

int DeltaKinematics_forward(DeltaKinematics *dk);
int DeltaKinematics_forward_angles(DeltaKinematics *dk, double thetaA, double thetaB, double thetaC);
int DeltaKinematics_inverse(DeltaKinematics *dk);
int DeltaKinematics_inverse_position(DeltaKinematics *dk, double x0, double y0, double z0);

int DeltaKinematics_delta_calcAngleYZ(double *Angle, double x0, double y0, double z0,
                                      double BassTri, double PlatformTri,
                                      double ArmLength, double RodLength);

int main() {
    DeltaKinematics DK;
    DeltaKinematics_setup(&DK, 50, 20, 75, 50);

    DK.x = 10;
    DK.y = 20;
    DK.z = 30;
    DeltaKinematics_inverse(&DK);
    printf("%.2f,%.2f,%.2f\n", DK.x, DK.y, DK.z);
    printf("%.2f,%.2f,%.2f\n", DK.a, DK.b, DK.c);

    return 0;
}

void DeltaKinematics_setup(DeltaKinematics *dk, double _ArmLength, double _RodLength, double _BassTri, double _PlatformTri) {
    dk->PlatformTri = _PlatformTri;
    dk->BassTri = _BassTri;
    dk->RodLength = _RodLength;
    dk->ArmLength = _ArmLength;
    dk->x = dk->y = dk->z = dk->a = dk->b = dk->c = 0.0;
}

int DeltaKinematics_forward(DeltaKinematics *dk) {
    return DeltaKinematics_forward_angles(dk, dk->a, dk->b, dk->c);
}

int DeltaKinematics_forward_angles(DeltaKinematics *dk, double thetaA, double thetaB, double thetaC) {
    dk->x = 0.0;
    dk->y = 0.0;
    dk->z = 0.0;
  
    double t = (dk->BassTri - dk->PlatformTri) * tan30 / 2.0;
    double dtr = pi / 180.0;

    thetaA *= dtr;
    thetaB *= dtr;
    thetaC *= dtr;

    double y1 = -(t + dk->ArmLength * cos(thetaA));
    double z1 = -dk->ArmLength * sin(thetaA);

    double y2 = (t + dk->ArmLength * cos(thetaB)) * sin30;
    double x2 = y2 * tan60;
    double z2 = -dk->ArmLength * sin(thetaB);

    double y3 = (t + dk->ArmLength * cos(thetaC)) * sin30;
    double x3 = -y3 * tan60;
    double z3 = -dk->ArmLength * sin(thetaC);

    double dnm = (y2 - y1) * x3 - (y3 - y1) * x2;

    double w1 = y1 * y1 + z1 * z1;
    double w2 = x2 * x2 + y2 * y2 + z2 * z2;
    double w3 = x3 * x3 + y3 * y3 + z3 * z3;

    double a1 = (z2 - z1) * (y3 - y1) - (z3 - z1) * (y2 - y1);
    double b1 = -((w2 - w1) * (y3 - y1) - (w3 - w1) * (y2 - y1)) / 2.0;

    double a2 = -(z2 - z1) * x3 + (z3 - z1) * x2;
    double b2 = ((w2 - w1) * x3 - (w3 - w1) * x2) / 2.0;

    double aV = a1 * a1 + a2 * a2 + dnm * dnm;
    double bV = 2.0 * (a1 * b1 + a2 * (b2 - y1 * dnm) - z1 * dnm * dnm);
    double cV = (b2 - y1 * dnm) * (b2 - y1 * dnm) + b1 * b1 + dnm * dnm * (z1 * z1 - dk->RodLength * dk->RodLength);

    double dV = bV * bV - 4.0 * aV * cV;
    if (dV < 0.0) {
        return non_existing_povar_error;
    }

    dk->z = -0.5 * (bV + sqrt(dV)) / aV;
    dk->x = (a1 * dk->z + b1) / dnm;
    dk->y = (a2 * dk->z + b2) / dnm;

    return no_error;
}

int DeltaKinematics_inverse(DeltaKinematics *dk) {
    return DeltaKinematics_inverse_position(dk, dk->x, dk->y, dk->z);
}

int DeltaKinematics_inverse_position(DeltaKinematics *dk, double x0, double y0, double z0) {
    dk->a = dk->b = dk->c = 0;
    int error = DeltaKinematics_delta_calcAngleYZ(&dk->a, x0, y0, z0, dk->BassTri, dk->PlatformTri, dk->ArmLength, dk->RodLength);
    if (error != no_error)
        return error;
    error = DeltaKinematics_delta_calcAngleYZ(&dk->b, x0 * cos120 + y0 * sin120, y0 * cos120 - x0 * sin120, z0,
                                              dk->BassTri, dk->PlatformTri, dk->ArmLength, dk->RodLength);
    if (error != no_error)
        return error;
    error = DeltaKinematics_delta_calcAngleYZ(&dk->c, x0 * cos120 - y0 * sin120, y0 * cos120 + x0 * sin120, z0,
                                              dk->BassTri, dk->PlatformTri, dk->ArmLength, dk->RodLength);
    return error;
}

int DeltaKinematics_delta_calcAngleYZ(double *Angle, double x0, double y0, double z0,
                                      double BassTri, double PlatformTri,
                                      double ArmLength, double RodLength) {
    double y1 = -0.5 * tan30 * BassTri;  // f/2 * tan(30 deg)
    y0 -= 0.5 * tan30 * PlatformTri;  // shift center to edge

    double aV = (x0 * x0 + y0 * y0 + z0 * z0 + ArmLength * ArmLength - RodLength * RodLength - y1 * y1) / (2.0 * z0);
    double bV = (y1 - y0) / z0;

    double dV = -(aV + bV * y1) * (aV + bV * y1) + ArmLength * (bV * bV * ArmLength + ArmLength);
    if (dV < 0) {
        return non_existing_povar_error;
    }

    double yj = (y1 - aV * bV - sqrt(dV)) / (bV * bV + 1);
    double zj = aV + bV * yj;
    *Angle = atan2(-zj, (y1 - yj)) * 180.0 / pi;

    return no_error;
}

