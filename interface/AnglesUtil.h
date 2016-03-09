#ifndef INC_ANGLESUTIL
#define INC_ANGLESUTIL
///////////////////////////////////////////////////////////////////////////////
// File: AnglesUtil.h
//
// Purpose: Provide useful functions for calculating angles and eta
//
// Created: 4-NOV-1998 Serban Protopopescu
// History: replaced old KinemUtil
// Modified: 23-July-2000 Add rapidity calculation
// 14-Aug-2003 converted all functions to double precision
///////////////////////////////////////////////////////////////////////////////
// Dependencies (#includes)

#include <iostream>
#include <math.h>
#include <stdlib.h> /* abs */

namespace kinem {
    const float ETA_LIMIT = 15.0;
    const float EPSILON = 1.E-10;
    
    /// Calculate phi from x, y
    inline double phi(double x, double y);
    /// Calculate phi for a line defined by xy1 and xy2 (xy2-xy1)
    inline double phi(double xy1[2], double xy2[2]);
    inline double phi(float xy1[2], float xy2[2]);
    
    // Calculate theta from x, y, z
    inline double theta(double x, double y, double z);
    // Calculate theta for a line defined by xyz1 and xyz2 (xyz2-xyz1)
    inline double theta(double xyz1[3], double xyz2[3]);
    inline double theta(float xyz1[3], float xyz2[3]);
    // Calculate theta from eta
    inline double theta(double etap);
    
    // Calculate eta from x, y, z (return also theta)
    inline double eta(double x, double y, double z);
    // Calculate eta for a line defined by xyz1 and xyz2 (xyz2-xyz1)
    inline double eta(double xyz1[3], double xyz2[3]);
    inline double eta(float xyz1[3], float xyz2[3]);
    // Calculate eta from theta
    inline double eta(double th);
    
    // Calculate rapidity from E, pz
    inline double y(double E, double pz);
    
    // Calculate phi1-phi2 keeping value between 0 and pi
    inline double delta_phi(double ph11, double phi2);
    // Calculate phi1-phi2 keeping value between -pi and pi
    inline double signed_delta_phi(double ph11, double phi2);
    // Calculate deltaR
    inline double delta_R(double eta1, double phi1, double eta2, double phi2);
    
    // Calculate unit std::vectors given two points
    inline void uvectors(double u[3], double xyz1[3], double xyz2[3]);
    inline void uvectors(float u[3], float xyz1[3], float xyz2[3]);
    
    inline double tanl_from_theta(double theta);
    inline double theta_from_tanl(double tanl);
}

inline double kinem::phi(double x, double y)
{
    double PHI = atan2(y, x);
    return (PHI >= 0.0) ? PHI : 2.0*M_PI + PHI;
}

inline double kinem::phi(double xy1[2], double xy2[2])
{
    double x = xy2[0] - xy1[0];
    double y = xy2[1] - xy1[1];
    return phi(x, y);
}

inline double kinem::phi(float xy1[2], float xy2[2])
{
    double dxy1[2] = { xy1[0], xy1[1] };
    double dxy2[2] = { xy2[0], xy2[1] };
    return phi(dxy1, dxy2);
}

inline double kinem::theta(double x, double y, double z)
{
    return atan2(sqrt(x*x + y*y), z);
}

inline double kinem::theta(double xyz1[3], double xyz2[3])
{
    double x = xyz2[0] - xyz1[0];
    double y = xyz2[1] - xyz1[1];
    double z = xyz2[2] - xyz1[2];
    return theta(x, y, z);
}

inline double kinem::theta(float xyz1[3], float xyz2[3])
{
    double dxyz1[3] = {xyz1[0], xyz1[1], xyz1[2]};
    double dxyz2[3] = {xyz2[0], xyz2[1], xyz2[2]};
    return theta(dxyz1, dxyz2);
}

inline double kinem::theta(double etap)
{
    return 2.0*atan(exp(-etap));
}

inline double kinem::eta(double x, double y, double z)
{
    double temp = sqrt(x*x + y*y + z*z) + EPSILON;
    return 0.5*log((temp + z) / (temp - z));
}

inline double kinem::eta(double xyz1[3], double xyz2[3])
{
    double x = xyz2[0] - xyz1[0];
    double y = xyz2[1] - xyz1[1];
    double z = xyz2[2] - xyz1[2];
    return eta(x, y, z);
}

inline double kinem::eta(float xyz1[3], float xyz2[3])
{
    double dxyz1[3] = {xyz1[0], xyz1[1], xyz1[2]};
    double dxyz2[3] = {xyz2[0], xyz2[1], xyz2[2]};
    return eta(dxyz1, dxyz2);
}

inline double kinem::eta(double th)
{
    if(th == 0) return ETA_LIMIT;
    if(th >= M_PI-0.0001) return -ETA_LIMIT;
    return -log(tan(th/2.0));
}

inline double kinem::y(double E, double pz)
{
    double temp = E + EPSILON;
    return 0.5 * log ((temp + pz)/(temp - pz));
}

inline double kinem::delta_phi(double phi1, double phi2)
{
    double PHI = abs(phi1 - phi2);
    return (PHI <= M_PI) ? PHI : 2.0*M_PI - PHI;
}

inline double kinem::signed_delta_phi(double phi1, double phi2)
{
    double phia = phi1;
    if(phi1 > M_PI) phia = phi1 - 2.0*M_PI;
    double phib = phi2;
    if(phi2 > M_PI) phib = phi2 - 2.0*M_PI;
    double dphi = phia - phib;
    if(dphi > M_PI) dphi -= 2.0*M_PI;
    if(dphi < -M_PI) dphi += 2.0*M_PI;
    return dphi;
}

inline double kinem::delta_R(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = kinem::delta_phi(phi1, phi2);
    return sqrt(deta*deta + dphi*dphi);
}

inline void kinem::uvectors(double u[3], double xyz1[3], double xyz2[3])
{
    double xdiff = xyz2[0] - xyz1[0];
    double ydiff = xyz2[1] - xyz1[1];
    double zdiff = xyz2[2] - xyz1[2];
    double s = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
    if(s > 0.0) {
        u[0]= xdiff/s;
        u[1]= ydiff/s;
        u[2]= zdiff/s;
    } else {
        u[0] = 0.0;
        u[1] = 0.0;
        u[2] = 0.0;
    }
}

inline void kinem::uvectors(float u[3], float xyz1[3], float xyz2[3])
{
    double du[3];
    double dxyz1[3] = {xyz1[0], xyz1[1], xyz1[2]};
    double dxyz2[3] = {xyz2[0], xyz2[1], xyz2[2]};
    uvectors(du, dxyz1, dxyz2);
    u[0] = du[0];
    u[1] = du[1];
    u[2] = du[2];
}

inline double kinem::tanl_from_theta(double theta)
{
    return tan(M_PI/2.0 - theta);
}

inline double kinem::theta_from_tanl(double tanl)
{
    return M_PI/2.0 - atan(tanl);
}

#endif // INC_ANGLESUTIL
