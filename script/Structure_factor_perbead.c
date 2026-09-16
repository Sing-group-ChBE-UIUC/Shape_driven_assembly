#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char** argv) {
    // ---- parameters ----
    const int warmup = 30;      // frames to skip
    const double qmin = 1e-3;
    const double qmax = 10;
    const int nmax = 30;         // integer range for 2π n / L
    const int cap_q = (nmax)*(nmax)*(nmax);
    const int nframe = 100;
    const int type_range = 3;

    // ---- files ----
    FILE *traj = fopen(argv[1], "r");
    char filename[64];
    // sprintf(filename, "S_q%d.txt", type_range);  // create the file name string
    sprintf(filename, "Sq_raw.txt");  // create the file name string
    FILE *skfile = fopen(filename, "w");          // open file for writing


    // ---- first frame header (orthogonal box assumed) ----
    int step, Ntot;
    double xlo,xhi,ylo,yhi,zlo,zhi, xko, yko, zko;
    char line[256];

    fscanf(traj,"ITEM: TIMESTEP\n%d\n", &step);
    fscanf(traj,"ITEM: NUMBER OF ATOMS\n%d\n", &Ntot);
    fgets(line, sizeof line, traj);                 // ITEM: BOX BOUNDS ...
    fscanf(traj,"%lf %lf\n%lf %lf\n%lf %lf\n", &xlo,&xhi, &ylo,&yhi, &zlo,&zhi);
    // fscanf(traj,"%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n", &xlo,&xhi,&xko, &ylo,&yhi,&yko, &zlo,&zhi, &zko);
    double Lx = xhi - xlo, Ly = yhi - ylo, Lz = zhi - zlo;
    fgets(line, sizeof line, traj);                 // ITEM: ATOMS ...

    // ---- buffers ----
    int *type = (int*)malloc(sizeof(int)*Ntot);
    double *x = (double*)malloc(sizeof(double)*Ntot);
    double *y = (double*)malloc(sizeof(double)*Ntot);
    double *z = (double*)malloc(sizeof(double)*Ntot);
    

    // first frame coords
    // for (int i=0, id, it; i<Ntot; i++) {
    //     fscanf(traj,"%d %d %lf %lf %lf\n", &id,&it,&x[i],&y[i],&z[i]);
    //     type[i] = it;
        // printf("%d\n", type[i]);
    for (int i=0, id, mol, it; i<Ntot; i++) {
        fscanf(traj,"%d %d %d %lf %lf %lf\n", &id,&mol, &it,&x[i],&y[i],&z[i]);
        // fscanf(traj,"%d %d %lf %lf %lf %d\n", &id,&it,&x[i],&y[i],&z[i],&mol);
        type[i] = it;
                // printf("%lf\n", x[i]);



    }
    // ---- build q-lattice and shell map (from first frame box) ----
    double *qx = (double*)malloc(sizeof(double)*cap_q);
    double *qy = (double*)malloc(sizeof(double)*cap_q);
    double *qz = (double*)malloc(sizeof(double)*cap_q);
    double    *qlist = (double*)malloc(sizeof(double)*cap_q);

    int nq = 0;
    for (int nx=0; nx<nmax; nx++) {
        for (int ny=0; ny<nmax; ny++) {
            for (int nz=0; nz<nmax; nz++) {
                if (nx==0 && ny==0 && nz==0) continue; // skip q=0
                double qxv = 2.0*M_PI*nx/Lx;
                double qyv = 2.0*M_PI*ny/Ly;
                double qzv = 2.0*M_PI*nz/Lz;
                double qabs = sqrt(qxv*qxv + qyv*qyv + qzv*qzv);
                if (qabs < qmin || qabs > qmax) continue;

                qx[nq] = qxv; qy[nq] = qyv; qz[nq] = qzv; qlist[nq] = qabs;
                nq++;
            }
        }
    }

    double *S_shell = (double*)calloc(cap_q, sizeof(double));
    long frames_used = 0;
    for(int frame_idx = 0; frame_idx < nframe; frame_idx++)
    {// ---- process frames (first frame is index 0) ----
 
        // count selected beads
        int Nsel = 0;
        // for (int i=0;i<Ntot;i++) if (type[i]==3 || type[i]==4) Nsel++; // Tapered
        for (int i=0;i<Ntot;i++) if (type[i] > 3) Nsel++; // Tapered
        // for (int i=0;i<Ntot;i++) if (type[i]>4 && type[i]<8) Nsel++; // double tapered
        if (frame_idx >= warmup)// && Nsel > 0) 
        {
            for (int s=1; s<cap_q; s++)
            {
                double re=0.0, im=0.0, qxv=qx[s], qyv=qy[s], qzv=qz[s];
                for (int i=0;i<Ntot;i++) 
                {
                    // if (type[i]==3 || type[i]==4) // Tapered
                    if (type[i]>3) // Tapered
                    // if (type[i]>4 && type[i]<8) // double tapered
                    {
                        double phase = qxv*x[i] + qyv*y[i] + qzv*z[i];
                        re += cos(phase);
                        im += sin(phase);

                        // printf("%d\n", type[i]);
                    }

                }
                double Sq = (re*re + im*im) / (double)Nsel;
                // double Sq = (re*re + im*im);
                S_shell[s] += Sq;
            } // Hi Monkey. this is jaccoon!
            
            frames_used++;
        }

        fscanf(traj,"ITEM: TIMESTEP\n%d\n", &step);
        fscanf(traj,"ITEM: NUMBER OF ATOMS\n%d\n", &Ntot);
        fgets(line, sizeof line, traj);                 // ITEM: BOX BOUNDS ...
        // fscanf(traj,"%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n", &xlo,&xhi,&xko, &ylo,&yhi,&yko, &zlo,&zhi, &zko);
        fscanf(traj,"%lf %lf\n%lf %lf\n%lf %lf\n", &xlo,&xhi,&ylo,&yhi,&zlo,&zhi);
        double Lx = xhi - xlo, Ly = yhi - ylo, Lz = zhi - zlo;
        fgets(line, sizeof line, traj);                 // ITEM: ATOMS ...
        // first frame coords
        for (int i=0, id, mol, it; i<Ntot; i++) {
            fscanf(traj,"%d %d %d %lf %lf %lf\n", &id,&mol, &it,&x[i],&y[i],&z[i]);
        // for (int i=0, id, it; i<Ntot; i++) {
        //     fscanf(traj,"%d %d %lf %lf %lf\n", &id, &it,&x[i],&y[i],&z[i]);
        //     type[i] = it;
        // for (int i=0, id, it, mol; i<Ntot; i++) {
        //     fscanf(traj,"%d %d %lf %lf %lf %d\n", &id,&it,&x[i],&y[i],&z[i],&mol);
                        type[i] = it;
                                        // printf("%lf\n", x[i]);

        }
    }


    for (int s=1; s<cap_q; s++) {
        // if (deg_shell[s] == 0 || frames_used == 0) continue;
        double qmid = qlist[s];
        double Savg = S_shell[s] / (double)frames_used;
        // printf("%lf %lf\n", qmid, Savg);
        fprintf(skfile, "%g %g\n", qmid, Savg);
    }

    fclose(traj);
    fclose(skfile);
    free(type); free(x); free(y); free(z);
    free(qx); free(qy); free(qz); free(qlist);
    free(S_shell); 
    return 0;
}
