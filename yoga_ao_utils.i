dtor = pi/180.;
radeg = 1./dtor;

RASC = 180*3600/pi;

func regressAx(y,x)
/* DOCUMENT coef = regressAx(y,x)
   returns the regression coefficient A, when fitting
   the data y by the function y=A.x
   This function expects that (y=Ax+B with B=0). This makes the linear fit more
   robust, since the regression line is constrained to go through (0,0).
 */
{
  return sum(y*x)/sum(x*x);
}

func fft_goodsize(size)
  /* DOCUMENT func fft_goodsize(size)
     find best size for a fft (radix 2, 3, 5 or 7) from size
   */
{
  aradix=[3,5,7];
  mradix = 2;

  for (i=1;i<=3;i++) {
    tmpf = log(size)/log(aradix(i));
    tmpl = long(tmpf);
    tmpf -= tmpl;
    mradix = ((tmpf > (log(size)/log(mradix) - long(log(size)/log(mradix))))(1) ? aradix(i) : mradix);
  }

  return mradix^(long(log(size)/log(mradix))+1);
}

func rotate(im,ang,cx,cy,zoom)
/* DOCUMENT imrot = rotate(im,ang,cx,cy,zoom)
            imrot = rotate(im,ang)

Rotates an image of an angle "ang" (in DEGREES).
The center of rotation is cx,cy.
A zoom factor can be applied.

(cx,cy) can be omitted :one will assume one rotates around the
center of the image.
If zoom is not specified, the default value of 1.0 is taken.

modif dg : allow to rotate a cube of images with one angle per image

 */
{
  if( is_void(zoom) ) zoom=1.0;
  if( zoom==0.0 ) zoom=1.0;
  if (numberof(ang) == 1)
    if(ang==0 & zoom==1.0) return im;
  
  ang *= pi/180.0;
  s = dimsof(im);
  nx = s(2);
  ny = s(3);
  if( is_void(cx) )
    cx = nx/2 + 1;
  if( is_void(cy) )
    cy = ny/2 + 1;
  x = indgen(1:nx)(,-:1:ny) - cx;
  y = indgen(1:ny)(-:1:nx,) - cy;
  x /= zoom;
  y /= zoom;
  
  if (numberof(ang)>1) {
    rxy = array(0.,nx,ny,2,numberof(ang));
    for (i=1;i<=numberof(ang);i++) {
      matrot = [[cos(ang(i)),-sin(ang(i))],[sin(ang(i)),cos(ang(i))]];
      rxy(,,,i) = [x,y](,,+) * matrot(,+) + [cx,cy](-,-,);
    }
  } else {
    matrot = [[cos(ang),-sin(ang)],[sin(ang),cos(ang)]];
    rxy = [x,y](,,+) * matrot(,+) + [cx,cy](-,-,);
  }
  
  nn = where(rxy<1);
  if( is_array(nn) )
    rxy(nn)=1.;
  
  if (numberof(ang) > 1) {
    rx = rxy(,,1,);
    ry = rxy(,,2,);
  } else {
    rx = rxy(,,1);
    ry = rxy(,,2);
  }
  nn = where(rx>(nx-1));
  if( is_array(nn) )
    rx(nn)=nx-1;
  nn = where(ry>(ny-1));
  if( is_array(nn) )
    ry(nn)=ny-1;

  wx = rx;
  wy = ry;
  rx = long(rx);   // partie entiere
  ry = long(ry);
  wx -= rx;        // partie fractionnaire
  wy -= ry;

  ind = rx + (ry-1)*nx;
  
  if (numberof(ang) > 1) {
    nim = indgen(numberof(ang))-1;
    nim *= (nx*ny);
    nim = nim(-:1:nx,-:1:ny,);
    ind += nim;
  }
  
  imr = (im(ind)*(1-wx) + im(ind+1)*wx)*(1-wy) + (im(ind+nx)*(1-wx)+im(ind+nx+1)*wx)*wy;
  return imr;
}


func fft_rotate(im,angle,xc=,yc=,gband=)
/* DOCUMENT fft_rotate(im,angle)
   im    : square image
   angle : rotation angle in degrees (clockwise)
   xc,yc : x,y positions of the rotation center
   
   high precision image rotation using fft
   no information loss if : image is shannon sampled
                            image has sufficiently large guard band
   using the algorithm from :
   "Fast Fourier method for the accurate rotation of sampled images"
   Kieran G. Larkin, Michael A. Oldfield, Hanno Klemm
   Optics Communications 139 (1997) 99-106

   routine by d. gratadour 31/05/2011
   SEE ALSO:
 */

{
  if (angle == 0.) return im;
  size = dimsof(im);
  if (size(2) != size(3)) error,"works only on square images";
  nx = size(2);

  if (angle >= 360) angle = angle % 360;
  if (angle > 180) return fft_rotate(im,angle-360,xc=xc,yc=yc,gband=gband);
  if (angle < -180) return fft_rotate(im,360+angle,xc=xc,yc=yc,gband=gband);

  if (gband != 0) {
    im2=array(double,2*nx,2*nx);
    if (nx %2 == 0)
      im2(nx-nx/2+1:nx+nx/2,nx-nx/2+1:nx+nx/2) = im;
    else
      im2(nx-nx/2:nx+nx/2,nx-nx/2:nx+nx/2) = im;      
    im = im2;
    nx *= 2;
  }
  
  if (xc == []) xc = ((nx/2)%2 == 0 ? nx/2+1 : nx/2);
  if (yc == []) yc = ((nx/2)%2 == 0 ? nx/2+1 : nx/2);

  theta = angle * pi/180.;
  
  if (angle > 90) theta = pi/2;
  if (angle < -90) theta = -pi/2;

  stepx = tan(theta/2);
  stepy = -1.*sin(theta);

  if ((nx/2)%2 == 0) {
    tiltx = -2.*pi/nx*(float(stepx)*(indgen(nx)-nx/2-0.5));
    tilty = -2.*pi/nx*(float(stepy)*(indgen(nx)-nx/2-0.5));
  } else {
    tiltx = -2.*pi/nx*(float(stepx)*(indgen(nx)-nx/2-1));
    tilty = -2.*pi/nx*(float(stepy)*(indgen(nx)-nx/2-1));
  }
  
  compltiltx=array(complex,nx);
  compltiltx.im=roll(tiltx);
  compltiltx = compltiltx(,-::nx-1);
  
  compltilty=array(complex,nx);
  compltilty.im=roll(tilty);
  compltilty = compltilty(-::nx-1,);
  
  col  = span(1,nx,nx)(,-:1:nx);
  lig  = transpose(col);
    
  tmpc=array(complex,nx,nx);

  tmpc = fft(exp(compltiltx*(lig-xc))*fft(im,[1,0]),[-1,0]);
  tmpc = fft(exp(compltilty*(col-yc))*fft(tmpc,[0,1]),[0,-1]);
  tmpc = fft(exp(compltiltx*(lig-xc))*fft(tmpc,[1,0]),[-1,0]);

  if (angle > 90) {
    return fft_rotate(tmpc.re/nx/nx/nx,angle-90.,gband=0);
  } else {
    if (angle < -90) {
      return fft_rotate(tmpc.re/nx/nx/nx,angle+90.,gband=0);
    } else {
      if ((nx/2)%2 == 0)
        return tmpc(nx/2-nx/4+1:nx/2+nx/4,nx/2-nx/4+1:nx/2+nx/4).re/nx/nx/nx;
      else
        return tmpc(nx/2-nx/4:nx/2+nx/4,nx/2-nx/4:nx/2+nx/4).re/nx/nx/nx;
    }
  }
}

func calc_dphi(phase,pup,den)
/* DOCUMENT
 * 
 */ 
{
  npix = dimsof(phase)(2);
  mi = p = dphi = array(float,[2,2*npix,2*npix]);
  mi(1:npix,1:npix) = phase;
  p(1:npix,1:npix)  = pup;

  //den  = (fft(fft(p,1)*conj(fft(p,1)),-1)).re;
  mask = den > max(den)*1.e-7;
  pix  = where(mask);
  //tmp = fft( abs(fft(mi*p))^2. ).re;
  //tmp(pix) /= (den)(pix);
  //dphi(pix) = tmp(1,1) - tmp(pix);
  dphi(pix) = fft(2*((fft(mi^2*p,1)*conj(fft(p,1))).re - abs(fft(mi*p,1))^2),-1).re(pix)/den(pix);

  return dphi; 
}

func calc_dphis(phase)
/* DOCUMENT
 * 
 */ 
{
  npix = dimsof(phase)(2);
  mi = p = dphi = array(float,[2,2*npix,2*npix]);
  mi(1:npix,1:npix) = phase;
  dphi = (mi - mi(1,1))^2;

  return dphi; 
}

func calc_dphif(phase,pup,den,conjftpup)
/* DOCUMENT 
 *  
 * KEYWORDS :  
 */ 
{
	npix = dimsof(phase)(2);
	mi = p = dphi = array(float,[2,2*npix,2*npix]);
	mi(1:npix,1:npix) = phase;
	p(1:npix,1:npix)  = pup;

	mask = den > max(den)*1.e-7;
	pix  = where(mask);
    //dphi(pix) = fft(2.*((fft(mi^2*p,1)*conjftpup) - abs(fft(mi*p,1))^2).re,-1).re(pix)/den(pix);
	dphi(pix) = fft(fft(mi^2*p, 1)*conjftpup + fft(p, 1)*conj(fft(mi^2*p, 1)) -2.*fft(mi*p, 1)*conj(fft(mi*p, 1)), -1).re(pix)/den(pix)
  
	return dphi;   
}

func circavg(a,center=,middle=)
/* DOCUMENT circavg
 *  
 * average=circavg(array[,center=,middle=])
 *
 * This routine returns the circular mean of an array. The routine dist is
 * used to compute an array of coordinate of points to be averaged.
 *
 * a      (input) : The array of which we want the circular mean. It can be
 *                  a long, a float or a double. Complex are not supported
 * center (input) : An array of the form [x,y] which give the position of
 *                  the origin for the circular mean calculation
 * middle (input) : A flag to indicate that the origin for the circular mean
 *                  calculation is the middle of the array a
 *
 * SEE ALSO: 
 */ 
{
  s=dimsof(a);

  if (!is_set(middle)) middle=0;
  if (s(1) != 2) write,"error - invalid dimensions";
  if (s(3) != s(2)) write,"error - invalid dimensions";

  dim=s(2);
  
  if (center!=[]) {
    s=dimsof(center);
    if ((s(1) != 1) | (s(1) != 2)) \
       write,"error - center has invalid dimensions";

    center=long(center);

    if (middle) {
      center=long([0,0]);
      write,"error - center and middle are not compatible keywords";
    }
  } else { 
    if (middle) center=long([0,0]);
    else center=[dim/2,dim/2];
  }
  
  r=long(roll(long(dist(dim)+.5)+1,[center(1),center(2)]));
  j=long(max(r));
  n=array(long,j);
  sx=array(double,j);
  dim2=long(dim)*long(dim);
  
  for (i=1;i<=dim2;i++) {
    j=r(i);
    sx(j)=sx(j)+a(i);
    n(j)=n(j)+1;
  }
  
  return sx/n;
}


func circavg_quad(a)
/* DOCUMENT circavg_quad
 *  
 * average=circavg_quad(array)
 *
 *
 * SEE ALSO: 
 */ 
{
  s=dimsof(a);

  if (s(1) != 2) write,"error - invalid dimensions";
  if (s(3) != s(2)) write,"error - invalid dimensions";

  dim=s(2);
  
  
  r=long(roll(dist(2*dim)+.5)+1)(1:dim,1:dim);
  j=long(max(r));
  n=array(long,j);
  sx=array(double,j);
  dim2=long(dim)*long(dim);
  
  for (i=1;i<=dim2;i++) {
    j=r(i);
    sx(j)=sx(j)+a(i);
    n(j)=n(j)+1;
  }
  
  return sx/n;
}

/*
 _____                     __   __ _    ___  
|  ___| __ ___  _ __ ___   \ \ / // \  / _ \ 
| |_ | '__/ _ \| '_ ` _ \   \ V // _ \| | | |
|  _|| | | (_) | | | | | |   | |/ ___ \ |_| |
|_|  |_|  \___/|_| |_| |_|   |_/_/   \_\___/ 
                                             
 */
func make_pupil(dim,pupd,xc=,yc=,real=,cobs=)
  /* DOCUMENT func make_pupil(dim,pupd,xc=,yc=,real=)
   */
{
  if (real == 1) {
    pup = exp(-(dist(dim,xc=xc,yc=yc)/(pupd/2.))^60.)^0.69314;
  } else {
    //    xc;yc;info,xc;
    //    tv,dist(dim,xc=xc,yc=yc);pause,2000;
    pup = dist(dim,xc=xc,yc=yc) < (pupd+1.)/2.;
  }
  if (is_set(cobs)) {
    if (real == 1) {
      pup -= exp(-(dist(dim,xc=xc,yc=yc)/(pupd*cobs/2.))^60.)^0.69314;
    } else {
      pup -= dist(dim,xc=xc,yc=yc) < (pupd*cobs+1.)/2.;
    }
  }
    
  return pup;
}
