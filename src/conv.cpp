// convert the data from "my" output to layout suitable for K.Werner
#include "conv.h"
#include "inc.h"

void convert_kw(const char *filename)
{
 double ****e, ****vx, ****vy, ****vz ;
 int nx, ny, nz, nt ;
 double xmin, xmax, ymin, ymax, zmin, zmax, tmin, tmax ;
 
 	ifstream finput(filename) ;
	if(!finput){
		cout << "error! ICs file is not opened\n" ;
		exit(1) ;
	}
	
	string dim_filename = filename ;
	dim_filename.append(".dim") ;
	ifstream dim_file (dim_filename.c_str()) ;

	// reading a header
	int nheadlines = 6 ;
	char ** header ;
	header = new char* [nheadlines] ;
	for(int i=0; i<nheadlines; i++){
		header[i] = new char [100] ;
		finput.getline(header[i], 100) ;
	}
	
	finput >> nt >> nx >> ny >> nz ;
	
	dim_file >> nt >> nx >> ny >> nz ;

	cout<<nt <<"  "<<nx<<"  "<<ny<<"  "<<nz<<endl ;
	e = new double*** [nt];
	vx = new double*** [nt];
	vy = new double*** [nt];
	vz = new double*** [nt];
	for(int it=0; it<nt; it++){
	e[it] = new double** [nx] ;
	vx[it] = new double** [nx] ;
	vy[it] = new double** [nx] ;
	vz[it] = new double** [nx] ;
	for(int ix=0; ix<nx; ix++){
		e[it][ix] = new double* [ny] ;
		vx[it][ix] = new double* [ny] ;
		vy[it][ix] = new double* [ny] ;
		vz[it][ix] = new double* [ny] ;
		for(int iy=0; iy<ny; iy++){
			e[it][ix][iy] = new double [nz] ;
			vx[it][ix][iy] = new double [nz] ;
			vy[it][ix][iy] = new double [nz] ;
			vz[it][ix][iy] = new double [nz] ;
		}
	}
	}
	double *t = new double [nt] ;
	double *x = new double [nx] ;
	double *y = new double [ny] ;
	double *z = new double [nz] ;
	
	for(int it=0; it<nt; it++){
	for(int ix=0; ix<nx; ix++)
	for(int iy=0; iy<ny; iy++)
	for(int iz=0; iz<nz; iz++){
		finput >> t[it] >> x[ix] >> y[iy] >> z[iz] >> e[it][ix][iy][iz] >> vx[it][ix][iy][iz] >> vy[it][ix][iy][iz] >> vz[it][ix][iy][iz]  ;
	}
	cout << t[it] << "         \r" << flush ;
	}
	cout << "done read\n" ;
	
	
//	for(int ix=0; ix<nx; ix++){
//	for(int iz=0; iz<nz; iz++){
//		cout<< x[ix] << "  " << vx[38][ix][ny/2][0] << "  " << e[38][ix][ny/2][0] << endl ;
//		cout<< z[iz] << "  " << vz[0][nx/2][ny/2][iz] << "  " << e[0][nx/2][ny/2][iz] << endl ;
//	}
	
/*	
 TCanvas *cc = new TCanvas("cc","cc");
 cc->GetPad(0)->SetFillColor(10) ;
 TGraph2D * g2d = new TGraph2D() ;
 int n=0; 
 for(int ix=0; ix<nx; ix++)
 //for(int iy=0; iy<ny; iy++){
 for(int iz=0; iz<nz; iz++){
//  g2d->SetPoint(n,ix,iy,v[14][ix][iy][0]) ;
	g2d->SetPoint(n,ix,iz,vz[0][ix][ny/2][iz]) ;
  n++ ;
  }
 
 g2d->GetXaxis()->SetTitle("x [fm]");
 g2d->GetXaxis()->CenterTitle();
 g2d->GetXaxis()->SetTitleOffset(1.8);
 g2d->GetYaxis()->SetTitle("y [fm]");
 g2d->GetYaxis()->CenterTitle();
 g2d->GetYaxis()->SetTitleOffset(1.5);
 g2d->GetZaxis()->SetTitle("\\epsilon [gev/fm3]");
 g2d->GetZaxis()->SetTitleOffset(1.2);
 gStyle->SetPalette(1) ;
 g2d->SetTitle("");
 g2d->Draw("surf1");
*/
 
 string filem = filename ;
 filem.append(".conv");
 ofstream fout(filem.c_str()) ;

	for(int i=0; i<nheadlines; i++){
		fout << header[i] << endl ;
	}
	
  fout<<nz<<"  "<<nt<<"  "<<nt<<"  "<<nx<<"  "<<ny<<endl ;
  fout<<z[0]<<"  "<<z[nz-1]<<"  "<<t[0]<<"  "<<t[nt-1]<<endl ;
  fout<<x[0]<<"  "<<x[nx-1]<<"  "<<y[0]<<"  "<<y[ny-1]<<endl ;
  
  for(int iy=0; iy<ny; iy++)
  for(int ix=0; ix<nx; ix++)
  for(int it=0; it<nt; it++)
  for(int iz=0; iz<nz; iz++)
	fout<< vx[it][ix][iy][iz] << "  " << vy[it][ix][iy][iz] << "  " << vz[it][ix][iy][iz] << "  " ;

// fout << endl ;
 for(int iy=0; iy<ny; iy++)
  for(int ix=0; ix<nx; ix++)
  for(int it=0; it<nt; it++)
  for(int iz=0; iz<nz; iz++)
	fout<< e[it][ix][iy][iz]<<"  " ;
	fout<<endl ;

}
