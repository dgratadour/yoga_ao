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
