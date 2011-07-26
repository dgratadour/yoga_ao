require,"yoga_aolib.i";
require,"util_fr.i";

func yoga_atmos_create(nscreen,r0,pupixsize,screen_size,frac,alt,windspeed,winddir,deltax,deltay,pupil)
/* DOCUMENT yoga_atmos_create
   g_atmos = yoga_atmos_create(nscreen,r0,pupixsize,screen_size,frac,alt,windspeed,winddir,deltax,deltay,pupil)

   creates an extrude ready yAtmos object on the gpu with all proper inits
   nscreen     : number of screens
   r0          : total r0 @ 0.5µm
   pupixsize   : pupil pixel size (in meters)
   screen_size : array of screen sizes (in pixels) for each layer
   frac     : array of r0 fractions per layers
   alt      : array of altitudes per layers
   wspeed   : array of wind speeds per layers
   wdir     : array of wind directions per layers
   deltax   : array of x displacement per iteration (one per layers)
   deltay   : array of y displacement per iteration (one per layers)
   pupil    : array containing the pupil

   SEE ALSO:
 */
{
  dirsave = YOGA_AO_SAVEPATH+"turbu/";
  mkdirp,dirsave;
  // creates dir if does not exists (else mkdirp does nothing)

 if (numberof(windspeed) != nscreen) error("wrong windspeed dimensions");
  if (numberof(winddir) != nscreen) error("wrong winddir dimensions");
  if (numberof(frac) != nscreen) error("wrong winddir dimensions");
  if (numberof(alt) != nscreen) error("wrong alt dimensions");
  //if (numberof(deltax) != nscreen) error("wrong deltax dimensions");
  //if (numberof(deltay) != nscreen) error("wrong deltay dimensions");
  if ((numberof(screen_size) != nscreen) && (numberof(screen_size) != 1)) {
    error("wrong screen_size dimensions");
  } else {
    if (numberof(screen_size) == 1)
      screen_size=screen_size(-::nscreen-1);
  }
  if (numberof(frac) == 1) frac = 1.0f;
  
  size2 = [];
  
  for (i=1;i<=nscreen;i++) {
    if ((!fileExist(swrite(format=dirsave+"A_%d.fits",screen_size(i)))) ||
        (!fileExist(swrite(format=dirsave+"B_%d.fits",screen_size(i)))) ||
        (!fileExist(swrite(format=dirsave+"istx_%d.fits",screen_size(i)))) ||
        (!fileExist(swrite(format=dirsave+"isty_%d.fits",screen_size(i))))) {
      AB, screen_size(i), A, B, ist;
      istx = ist;
      test = array(0,screen_size(i),screen_size(i));
      test(istx) = indgen(dimsof(A)(3));
      test2=transpose(test);
      isty = sort(test2(*))(screen_size(i)*screen_size(i)-dimsof(A)(3)+1:);
      fits_write,swrite(format=dirsave+"A_%d.fits",screen_size(i)),A;
      fits_write,swrite(format=dirsave+"B_%d.fits",screen_size(i)),B;
      fits_write,swrite(format=dirsave+"istx_%d.fits",screen_size(i)),istx;
      fits_write,swrite(format=dirsave+"isty_%d.fits",screen_size(i)),isty;
    } else {
      A = fits_read(swrite(format=dirsave+"A_%d.fits",screen_size(i)));
      B = fits_read(swrite(format=dirsave+"B_%d.fits",screen_size(i)));
      istx = fits_read(swrite(format=dirsave+"istx_%d.fits",screen_size(i)));
      isty = fits_read(swrite(format=dirsave+"isty_%d.fits",screen_size(i)));
    }
    grow,size2,dimsof(A)(3);
  }
  
  /*
    f(i) = (r0tot/r0(i))^(5./3)
    r0(i) = r0tot / f(i)^(3./5)
    here r0tot must be given in pixels @ 0.5µm
  */

  r0 = r0  / frac^(3./5.);
  atmos_obj = yoga_atmos(int(nscreen),float(r0/pupixsize),long(screen_size),long(size2),float(alt),float(windspeed),float(winddir),float(deltax),float(deltay),float(pupil));
  
  atmos_obj;

  for (i=1;i<=nscreen;i++) {
    A = fits_read(swrite(format=dirsave+"A_%d.fits",screen_size(i)));
    B = fits_read(swrite(format=dirsave+"B_%d.fits",screen_size(i)));
    istx = fits_read(swrite(format=dirsave+"istx_%d.fits",screen_size(i)));
    isty = fits_read(swrite(format=dirsave+"isty_%d.fits",screen_size(i)));
    init_tscreen,atmos_obj,alt(i),float(A),float(B),int(istx-1),int(isty-1),int(1234*i);
    tic;
    for (cc=1;cc<=2*screen_size(i);cc++) extrude_tscreen,atmos_obj,alt(i);
    write,format="extrude time : %f\n",tac()/(2*screen_size(i));
  }
  
  return atmos_obj;
}


func create_screen(r0,pupixsize,screen_size,&A,&B,&ist)
/* DOCUMENT create_screen
   screen = create_screen(r0,pupixsize,screen_size,&A,&B,&ist)

   creates a phase screen and fill it with turbulence
   r0          : total r0 @ 0.5µm
   pupixsize   : pupil pixel size (in meters)
   screen_size : screen size (in pixels)
   A           : A array for future extrude
   B           : B array for future extrude
   ist         : istencil array for future extrude

   SEE ALSO:
 */
{
  AB, screen_size, A, B, ist;   // initialisation for A and B matrices for phase extrusion
  screen = array(0.0,screen_size,screen_size);   // init of first phase screen
  for(i=1;i<=2*screen_size;i++) screen=extrude(screen, r0/pupixsize, A, B, ist);
  return screen;
}

/*  *************************************************************
Extrude functions.
Initially written by e. gendron
modified by d. gratadour for compatibility with large screen sizes
************************************************************* */


func createStencil(n, &Z_x, &Z_y, &X_x, &X_y, &istencil)
/* DOCUMENT createStencil, n, Z_x, Z_y, X_x, X_y, istencil
     
   SEE ALSO:
 */
{ 
  Z_x = indgen(1:n)(,-:1:n);
  Z_y = indgen(1:n)(-:1:n,);
  
  X_x = array(n+1,n);
  X_y = indgen(1:n);
  
  // creation stencil
  ns = long( (log((n+1))/log(2))+1 );
  stencil = array(0,n,n);
  stencil(1,)=1;
  for(i=2;i<=ns+1;i++) {
    //stencil(i,::2^(i-1)) = 1;
    stencil(2^(i-2)+1,::2^(i-1)) = 1;
    stencil(2^(i-1),1) = 1;
  }
  stencil = roll(stencil,[0,n/2]);
  stencil = stencil(::-1,);

  istencil = where(stencil);
}





func Cxz(n, Z_x, Z_y, X_x, X_y, istencil)
/* DOCUMENT xz = Cxz(n, istencil)
   Cxz computes the covariance matrix between the new phase vector x (new
   column for the phase screen), and the already known phase values z.

   The known values z are the values of the phase screen that are pointed by
   the stencil indexes (istencil)

   SEE ALSO: Czz Cxx
 */
{
  xz = -((X_x(,-)-Z_x(istencil)(-,))^2 + (X_y(,-)-Z_y(istencil)(-,))^2)^(5./6.) +
    (((Z_x(n)-X_x)^2 + (Z_y(n)-X_y)^2)^(5./6.))(,-) +
    (((Z_x(n)-Z_x(istencil))^2 + (Z_y(n)-Z_y(istencil))^2)^(5./6.))(-,);
  xz *= 0.5 * 6.88;
  return xz;
}





func Cxx(n, X_x, X_y)
/* DOCUMENT xx = Cxx(X_x, X_y)
   Cxx computes the covariance matrix of the new phase vector x (new
   column for the phase screen).


   SEE ALSO: Czz Cxz
 */
{
  
  xx = -((X_x(,-)-X_x(-,))^2 + (X_y(,-)-X_y(-,))^2)^(5./6.) +
    (((Z_x(n)-X_x)^2 + (Z_y(n)-X_y)^2)^(5./6.))(,-) +
    (((Z_x(n)-X_x)^2 + (Z_y(n)-X_y)^2)^(5./6.))(-,);
  xx *= 0.5 * 6.88;
  return xx;
}





func Czz(n, Z_x, Z_y, istencil)
/* DOCUMENT zz = Czz(n, Z_x, Z_y, istencil)
   Czz computes the covariance matrix of the already known phase values z.

   The known values z are the values of the phase screen that are pointed by
   the stencil indexes (istencil)

   SEE ALSO: AB createStencil Cxz Cxx extrude
 */
{  
  zz = -((Z_x(istencil)(,-)-Z_x(istencil)(-,))^2 + (Z_y(istencil)(,-)-Z_y(istencil)(-,))^2)^(5./6.) +
    (((Z_x(n)-Z_x(istencil))^2 + (Z_y(n)-Z_y(istencil))^2)^(5./6.))(,-) +
    (((Z_x(n)-Z_x(istencil))^2 + (Z_y(n)-Z_y(istencil))^2)^(5./6.))(-,);
  zz *= 0.5 * 6.88;
  return zz;
}





func AB(n, &A, &B, &istencil)
/* DOCUMENT AB, n, A, B, istencil
     This function initializes some matrices A, B and a list of stencil indexes
     istencil for iterative extrusion of a phase screen.

     The method used is described by Fried & Clark in JOSA A, vol 25, no 2, p463, Feb 2008.
     The iteration is :
     x = A(z-zRef) + B.noise + zRef
     with z a vector containing "old" phase values from the initial screen, that are listed
     thanks to the indexes in istencil.
     
   SEE ALSO: extrude createStencil Cxx Cxz Czz
 */
{
  createStencil, n, Z_x, Z_y, X_x, X_y, istencil;   // init of matrices A and B, and istencil
  zz = Czz(n, Z_x, Z_y, istencil);                  // compute cov matrices
  xz = Cxz(n, Z_x, Z_y, X_x, X_y, istencil);
  xx = Cxx(n, X_x, X_y);

  //  zz_1 = LUsolve( zz );    // zz cannont be inverted because of the reference zRef, which makes it singular.
  s = s1 = SVdec(zz,uuu,vt);   // SV decomp of zz
  s1(0)=1;
  s1 = 1./s1;
  s1(0)=0;                     // the null eignevalue is not inverted
  zz_1 =  (uuu*s1(-,))(,+) * vt(+,);   // inversion, with the null eigenvalue left to 0
  
  A = xz(,+) * zz_1(+,);
  
  bbt = xx - A(,+)*xz(,+);
  l = SVdec(bbt,uu);
  B = uu*sqrt(l)(-,);
}




func extrude(p,r0,A,B,istencil)
/* DOCUMENT p1 = extrude(p,r0,A,B,istencil)
     Extrudes a phase screen p1 from initial phase screen p.
     p1 prolongates p by 1 column on the right end.
     r0 is expressed in pixels
     
     The method used is described by Fried & Clark in JOSA A, vol 25, no 2, p463, Feb 2008.
     The iteration is :
     x = A(z-zRef) + B.noise + zRef
     with z a vector containing "old" phase values from the initial screen, that are listed
     thanks to the indexes in istencil.

     Examples
     n = 32;
     AB, n, A, B, istencil;
     p = array(0.0,n,n);
     p1 = extrude(p,r0,A,B,istencil);
     pli, p1
     
   SEE ALSO: AB() createStencil() Cxx() Cxz() Czz()
 */
{
  amplitude = r0^(-5./6);
  n = dimsof(p)(2);
  z = p(istencil);
  zref = p(n);
  z -= zref;
  newColumn = A(,+) * z(+) + B(,+)*(random_n(n)* amplitude )(+) + zref;
  p1 = array(0.0,n,n);
  p1(1:n-1,)=p(2:n,);
  p1(n,) = newColumn;
  return p1;
}







/*
  ___  _     _      _    ____ ___ 
 / _ \| | __| |    / \  |  _ \_ _|
| | | | |/ _` |   / _ \ | |_) | | 
| |_| | | (_| |  / ___ \|  __/| | 
 \___/|_|\__,_| /_/   \_\_|  |___|
                                  
*/
/*
func yoga_extrude_init(screen_size,&d_A,&d_B,&d_screen,&d_z,&d_ytmp,&d_noise,&d_tmpscr,A=,B=,istencil=)
{
  //screen_size = 128;
  "Generating Matrices ..";
  if ((A == []) || (B == []) || (istencil == []))
    AB, screen_size, A, B, istencil;
  // initialisation for A and B matrices for phase extrusion
  "Allocating GPU mem ..";
  d_A = yoga_obj(float(A));
  d_B = yoga_obj(float(B));
  d_screen = yoga_obj(array(0.0f,screen_size,screen_size));
  d_tmpscr = yoga_obj(array(0.0f,screen_size-1,screen_size));
  d_z = yoga_compact(d_screen,int(istencil-1));
  d_noise = yoga_random_n("float",[1,screen_size]);
  d_ytmp = yoga_obj("float",[1,screen_size]);
}

func yoga_extrude(&d_screen,r0,d_A,d_B,d_z,d_tmpscr,d_noise,d_ytmp,screen_size)
{
  profile = [];
  amplitude = r0^(-5./6);
  stime = tic();
  yoga_compact,d_screen,d_z;
  grow,profile,tac(stime);
  //zref = yoga_getvalue(d_screen,screen_size);
  //zref=0.0f;
  //grow,profile,tac(stime);
  yoga_plusai,d_z,d_screen,int(screen_size),int(-1);
  grow,profile,tac(stime);
  yoga_mv,d_ytmp,d_A,d_z;
  grow,profile,tac(stime);
  yoga_random_n,d_noise;
  grow,profile,tac(stime);
  yoga_mv,d_ytmp,d_B,d_noise,amplitude,1.0;
  grow,profile,tac(stime);
  yoga_plusai,d_ytmp,d_screen,int(screen_size);
  grow,profile,tac(stime);
  yoga_getarray,d_tmpscr,d_screen,2:screen_size,;
  grow,profile,tac(stime);
  yoga_fillarray,d_screen,d_tmpscr,1:screen_size-1,; //can use a memcopy
  grow,profile,tac(stime);
  yoga_fillarray,d_screen,d_ytmp,screen_size:screen_size,;
  grow,profile,tac(stime);
  return profile;
}

func yoga_extrude_fast(&d_screen,r0,d_A,d_B,d_z,d_tmpscr,d_noise,d_ytmp,d_anoise,screen_size,iter)
{
  profile = [];
  amplitude = r0^(-5./6);
  stime = tic();
  yoga_compact,d_screen,d_z;
  grow,profile,tac(stime);
  //zref = yoga_getvalue(d_screen,screen_size);
  //zref=0.0f;
  //grow,profile,tac(stime);
  yoga_plusai,d_z,d_screen,int(screen_size),int(-1);
  grow,profile,tac(stime);
  yoga_mv,d_ytmp,d_A,d_z;
  grow,profile,tac(stime);
  //yoga_random_n,d_noise;
  yoga_getarray,d_noise,d_anoise,1:screen_size,iter:iter;
  grow,profile,tac(stime);
  yoga_mv,d_ytmp,d_B,d_noise,amplitude,1.0;
  grow,profile,tac(stime);
  yoga_plusai,d_ytmp,d_screen,int(screen_size);
  grow,profile,tac(stime);
  yoga_getarray,d_tmpscr,d_screen,2:screen_size,;
  grow,profile,tac(stime);
  yoga_fillarray,d_screen,d_tmpscr,1:screen_size-1,; //can use a memcopy
  grow,profile,tac(stime);
  yoga_fillarray,d_screen,d_ytmp,screen_size:screen_size,;
  grow,profile,tac(stime);
  return profile;
}

func yoga_create_screen(r0,pupixsize,screen_size)
{
  
  yoga_extrude_init,screen_size,d_A,d_B,d_screen,d_z,d_ytmp,d_noise,d_tmpscr;

  profile = 0.0f;
  for(i=1;i<=2*screen_size;i++)
    profile += yoga_extrude(d_screen,r0/pupixsize,d_A,d_B,d_z,d_tmpscr,d_noise,d_ytmp,screen_size);
  profile /= (2*screen_size);
  return screen;
}

func turb_bench(r0,pupixsize,screen_size)
{
  tic;
  "CPU inits";
  if ((!fileExist(swrite(format="A_%d.fits",screen_size))) ||
      (!fileExist(swrite(format="B_%d.fits",screen_size))) ||
      (!fileExist(swrite(format="ist_%d.fits",screen_size))))
    AB, screen_size, A, B, ist;
  else {
    A = fits_read(swrite(format="A_%d.fits",screen_size));
    B = fits_read(swrite(format="B_%d.fits",screen_size));
    ist = fits_read(swrite(format="ist_%d.fits",screen_size));
  }
  
  if (!fileExist(swrite(format="A_%d.fits",screen_size)))
      fits_write,swrite(format="A_%d.fits",screen_size),A;
  if (!fileExist(swrite(format="B_%d.fits",screen_size)))
    fits_write,swrite(format="B_%d.fits",screen_size),B;
  if (!fileExist(swrite(format="ist_%d.fits",screen_size)))
    fits_write,swrite(format="ist_%d.fits",screen_size),ist;
  
  screen1 = array(0.0,screen_size,screen_size);   // init of first phase screen
  tac();
  
  "GPU inits";
  tic;
  yoga_extrude_init,screen_size,d_A,d_B,d_screen,d_z,d_ytmp,d_noise,d_tmpscr,A=A,B=B,istencil=ist;
  tac();

  profile = 0.0f;
  for(i=1;i<=2*screen_size;i++) {
    write,format="\r iter # %d",i;
    profile += yoga_extrude(d_screen,r0/pupixsize,d_A,d_B,d_z,d_tmpscr,d_noise,d_ytmp,screen_size);
  }
  write,format="\ngpu time : %f\n",profile(0)/(2*screen_size);
  
  tic;
  for(i=1;i<=100;i++) {
    write,format="\r iter # %d",i;
    screen1 = extrude(screen1, r0/pupixsize, A, B, ist);
  }
  write,format="\ncpu time : %f\n",tac()/(100);

  screen2 = d_screen();
  screen1 -= avg(screen1);
  screen2 -= avg(screen2);
  

  write,"gpu screen : min max avg rms : %f %f %f %f\n",min(screen2(*)),max(screen2(*)),screen2(*)(avg),screen2(*)(rms);
  write,"cpu screen : min max avg rms : %f %f %f %f\n",min(screen1(*)),max(screen1(*)),screen1(*)(avg),screen1(*)(rms);
  error;
}

func turb_bench_fast(r0,pupixsize,screen_size)
{
  tic;
  "CPU inits";
  if ((!fileExist(swrite(format="A_%d.fits",screen_size))) ||
      (!fileExist(swrite(format="B_%d.fits",screen_size))) ||
      (!fileExist(swrite(format="ist_%d.fits",screen_size))))
    AB, screen_size, A, B, ist;
  else {
    A = fits_read(swrite(format="A_%d.fits",screen_size));
    B = fits_read(swrite(format="B_%d.fits",screen_size));
    ist = fits_read(swrite(format="ist_%d.fits",screen_size));
  }
  
  if (!fileExist(swrite(format="A_%d.fits",screen_size)))
      fits_write,swrite(format="A_%d.fits",screen_size),A;
  if (!fileExist(swrite(format="B_%d.fits",screen_size)))
    fits_write,swrite(format="B_%d.fits",screen_size),B;
  if (!fileExist(swrite(format="ist_%d.fits",screen_size)))
    fits_write,swrite(format="ist_%d.fits",screen_size),ist;
  
  screen1 = array(0.0,screen_size,screen_size);   // init of first phase screen
  tac();
  
  "GPU inits";
  tic;
  yoga_extrude_init,screen_size,d_A,d_B,d_screen,d_z,d_ytmp,d_noise,d_tmpscr,A=A,B=B,istencil=ist;
  tac();

  "Generate noise on GPU";
  tic;
  d_anoise = yoga_random_n("float",[2,screen_size,2*screen_size]);
  tac();

  profile = 0.0f;
  for(i=1;i<=2*screen_size;i++) {
    write,format="\r iter # %d",i;
    profile += yoga_extrude_fast(d_screen,r0/pupixsize,d_A,d_B,d_z,d_tmpscr,d_noise,d_ytmp,d_anoise,screen_size,i);
  }
  write,format="\ngpu time : %f\n",profile(0)/(2*screen_size);
  
  tic;
  for(i=1;i<=2*screen_size;i++) {
    write,format="\r iter # %d",i;
    screen1 = extrude(screen1, r0/pupixsize, A, B, ist);
  }
  write,format="\ncpu time : %f\n",tac()/(2*screen_size);

  screen2 = d_screen();
  screen1 -= avg(screen1);
  screen2 -= avg(screen2);
  

  write,"gpu screen : min max avg rms : %f %f %f %f\n",min(screen2(*)),max(screen2(*)),screen2(*)(avg),screen2(*)(rms);
  write,"cpu screen : min max avg rms : %f %f %f %f\n",min(screen1(*)),max(screen1(*)),screen1(*)(avg),screen1(*)(rms);
  error;
}

r0 = float(0.16);pupixsize = float(0.03125);screen_size=512;
turb_bench,r0,pupixsize,screen_size;
*/
/*
Profiles on a c2050 :
curand host api
4096 : accel factor ~30
  compact  /  -zref    /    mv1    /   rng    /  mv2    /  +zref   /  get  /   fill  /   fill
[7.31816e-05,0.000120031,0.000134508,0.0146606,0.0147031,0.0147558,0.0148052,0.0148527,0.0148986]

          [4.68496e-05,1.44769e-05,0.0145261,4.24513e-05,5.27482e-05,4.93287e-05,4.75156e-05,4.59278e-05]

1024 : accel factor ~15
[2.81489e-05,5.71714e-05,6.52273e-05,0.00127357,0.00128396,0.00131188,0.00133945,0.00136672,0.00139357]

          [2.90226e-05,8.05582e-06,0.00120834,1.03838e-05,2.79212e-05,2.75772e-05,2.72667e-05,2.68492e-05]

curand device api
1024 : accel factor ~15
[4.72562e-05,9.3541e-05,0.000107719,0.00122242,0.00124069,0.00128585,0.00133204,0.00137791,0.00142375]

           [4.62848e-05,1.41779e-05,0.0011147,1.82661e-05,4.5159e-05,4.61912e-05,4.58696e-05,4.58435e-05]

4096 : accel factor ~30
[4.91733e-05,7.97563e-05,8.85014e-05,0.0147314,0.0147655,0.0148007,0.0148314,0.0148608,0.0148887]

           [3.05829e-05,8.74512e-06,0.0146429,3.40792e-05,3.52175e-05,3.07083e-05,2.94099e-05,2.78721e-05]


conclusion : our kernel launch does a great job !
better to stick to host api while we can

*/
