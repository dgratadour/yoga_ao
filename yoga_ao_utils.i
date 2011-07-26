dtor = pi/180.;
radeg = 1./dtor;

RASC = 180*3600/pi;

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
    im2(nx-nx/2:nx+nx/2,nx-nx/2:nx+nx/2) = im;
    im = im2;
    nx *= 2;
  }
  
  if (xc == []) xc = ((nx/2)%2 == 0 ? nx/2+0.5 : nx/2);
  if (yc == []) yc = ((nx/2)%2 == 0 ? nx/2+0.5 : nx/2-1);

  theta = angle * pi/180.;
  
  if (angle > 90) theta = pi/2;
  if (angle < -90) theta = -pi/2;

  stepx = tan(theta/2);
  stepy = -1.*sin(theta);

  // if (nx/2%2) {
  tiltx = -2.*pi/nx*(float(stepx)*(indgen(nx)-nx/2-0.5));
  tilty = -2.*pi/nx*(float(stepy)*(indgen(nx)-nx/2-0.5));
    //  } else {
  //tiltx = -2.*pi/nx*(float(stepx)*(indgen(nx)-nx/2-1));
  //tilty = -2.*pi/nx*(float(stepy)*(indgen(nx)-nx/2-1));
    //  }
  
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
  dphi = (mi - mi(1,1))^2;
  /*
  p(1:npix,1:npix)  = pup;

  //den  = (fft(fft(p,1)*conj(fft(p,1)),-1)).re;
  mask = den > max(den)*1.e-7;
  pix  = where(mask);
  tmp = fft( abs(fft(mi*p))^2. ).re;
  tmp(pix) /= (den)(pix);
  dphi(pix) = tmp(1,1) - tmp(pix);
  */
  //dphi(pix) = fft(2*((fft(mi^2*p,1)*conj(fft(p,1))).re - abs(fft(mi*p,1))^2),-1).re(pix)/den(pix);

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
