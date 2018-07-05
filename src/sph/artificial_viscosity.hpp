#ifndef ARTIFICIAL_VISCOSITY_HPP
#define ARTIFICIAL_VISCOSITY_HPP
using namespace std;

void art_force(int ntotal,int niac,vector<int> pair_i,vector<int> pair_j,
               vector <double> position,vector <double> velocity, vector <double> &art_acceleration, vector <double> &art_energy,
               double hsml, double mass, vector <double> rho,
               vector <double> pressure,
               vector <double> dwdx)
{
    int i,j,k;
    double xij,vx_ij,xpsv,xpsx,muv,piv,hx, mrho, mc, mhsml;
    double qa = 1.0;
    double qb = 1.0;
    double etq = 1.0e-1;
    for(k=0;k<ntotal;k++)
    {
        art_acceleration[k] = 0.0;
        art_energy[k] = 0.0;
    }
    for(k=0;k<niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];

        xij = position[i] - position[j];
        vx_ij = velocity[i] - velocity[j];

        xpsv = vx_ij * xij;
        xpsx = xij * xij;

        if(xpsv<0.0)
        {
            mhsml = hsml;//0.5 * (hsml[i] + hsml[j]);
            muv = mhsml * xpsv/(xpsx + mhsml * mhsml * etq * etq);
            mc = (sqrt(1.4 * pressure[i]/rho[i]) + sqrt(1.4 * pressure[j]/rho[j])) * 0.5;
            mrho = (rho[i] + rho[j]) * 0.5;
            piv = (qb * muv - qa * mc) * muv/mrho;
            hx = - piv * dwdx[k];

            art_acceleration[i] = art_acceleration[i] + mass * hx;
            art_acceleration[j] = art_acceleration[j] - mass * hx;

            art_energy[i] = art_energy[i] + 0.5 * (velocity[i] - velocity[j]) * hx * mass;
            art_energy[j] = art_energy[j] + 0.5 * (velocity[i] - velocity[j]) * hx * mass;
        }
    }

}

#endif // ARTIFICIAL_VISCOSITY_HPP

