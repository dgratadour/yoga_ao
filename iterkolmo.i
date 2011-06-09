
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





//  *************************************************************



func createStencilElongated(nx, ny, elong, &Z_x, &Z_y, &X_x, &X_y, &istencil)
/* DOCUMENT createStencil, n, Z_x, Z_y, X_x, X_y, istencil
     
   SEE ALSO:
 */
{
  elong = long(elong);
  Z_x = indgen(1:nx)(,-:1:ny) / double(elong);
  Z_y = indgen(1:ny)(-:1:nx,);
  
  X_x = array(nx+1,ny) / double(elong);
  X_y = indgen(1:ny);
  
  // creation stencil
  ns = long( (log((ny+1))/log(2))+1 );
  stencil = array(0,nx,ny);
  stencil(1:2,)=1;
  /*
  for(i=2;i<=ns+1;i++) {
    // stencil(i,::2^(i-1)) = 1;   // formule originale sans elongation
    // stencil((i-2)*elong+2,::2^(i-1)) = 1;   // essai
    //stencil(i*(i-1)+1,::2^(i-1)) = 1;  // on voit toujours la marque du stencil...
    stencil(i*(i-1)+1,) = 1;
  }
  */
  // autre méthode
  stencil(1:2,)=1;
  i=2;
  while( (i*(i-1)+1)<nx) {
    stencil(i*(i-1)+1,) = 1;
    i++;
  }
  
  stencil = roll(stencil,[0,ny/2]);
  stencil = stencil(::-1,);
  
  istencil = where(stencil);
}




func ABElongated(nx, ny, elong, &A, &B, &istencil)
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
  createStencilElongated, nx, ny, elong, Z_x, Z_y, X_x, X_y, istencil;   // init of matrices A and B, and istencil
  zz = Czz(nx, Z_x, Z_y, istencil);                  // compute cov matrices
  xz = Cxz(nx, Z_x, Z_y, X_x, X_y, istencil);
  xx = Cxx(nx, X_x, X_y);

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





func extrudeElongated(p,r0,A,B,istencil)
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
  nx = dimsof(p)(2);
  ny = dimsof(p)(3);
  z = p(istencil);
  zref = p(nx);
  z -= zref;
  newColumn = A(,+) * z(+) + B(,+)*(random_n(ny)* amplitude )(+) + zref;
  p1 = array(0.0,nx,ny);
  p1(1:nx-1,)=p(2:nx,);
  p1(nx,) = newColumn;
  return p1;
}



