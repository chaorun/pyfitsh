/*****************************************************************************/
/* gropt_interface.c                                                         */
/* C interface for gropt geometrical optics (no FILE I/O)                    */
/* Extracted from gropt.c, all fprintf replaced with strbuf or arrays        */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "optcalc.h"
#include "gropt_interface.h"
#include "list.h"
#include "math/polygon.h"

/*****************************************************************************/
/* strbuf: dynamic string buffer (replaces FILE* fprintf)                    */
/*****************************************************************************/

typedef struct
 {	char	*data;
	size_t	len;
	size_t	cap;
 } gropt_strbuf;

static void gropt_strbuf_init(gropt_strbuf *sb)
{
 sb->data=NULL;
 sb->len=0;
 sb->cap=0;
}

static void gropt_strbuf_appendf(gropt_strbuf *sb,const char *fmt,...)
{
 va_list	ap;
 int		needed;

 va_start(ap,fmt);
 needed=vsnprintf(NULL,0,fmt,ap);
 va_end(ap);

 if ( sb->len+needed+1 > sb->cap )
  {	sb->cap=(sb->len+needed+1)*2;
	if ( sb->cap<1024 )	sb->cap=1024;
	sb->data=(char *)realloc(sb->data,sb->cap);
  }

 va_start(ap,fmt);
 vsnprintf(sb->data+sb->len,needed+1,fmt,ap);
 va_end(ap);

 sb->len+=needed;
}

/*****************************************************************************/
/* gropt_transfer: ray transfer matrix analysis                              */
/*****************************************************************************/

int gropt_transfer(struct optics *opt,double wavelength,
	double *out_focal_plane,double *out_eff_focus)
{
 transfer_matrix	m;
 double			zend,distance,eff_focus,z_focus;

 optcalc_compute_transfer_matrix(opt,wavelength,m,0.0,&zend);
 distance=-m[0][0]/m[1][0];
 eff_focus=m[0][1]+distance*m[1][1];
 z_focus=zend+distance;

 if ( out_focal_plane != NULL )	*out_focal_plane=z_focus;
 if ( out_eff_focus != NULL )	*out_eff_focus=eff_focus;

 return(0);
}

/*****************************************************************************/
/* gropt_spot_diagram: spot diagram analysis                                 */
/*****************************************************************************/

int gropt_spot_diagram(struct optics *opt,
	double wavelength,double aperture_radius,int nrings,
	double angle_nx,double angle_ny,double zstart,double pixel_scale,
	double *out_xy,int out_max,int *out_n,
	double *out_center_x,double *out_center_y)
{
 int	i,j,nsec,cnt;
 double	sx0,sy0;
 vector	v0,n0;
 double	angle_nz;

 angle_nz=sqrt(1.0-(angle_nx*angle_nx+angle_ny*angle_ny));

 optcalc_glass_refraction_precompute(opt,wavelength);

 v0[0]=0.0;
 v0[1]=0.0;
 v0[2]=zstart-0.0001;
 n0[0]=angle_nx;
 n0[1]=angle_ny;
 n0[2]=angle_nz;
 optcalc_ray_trace(opt,wavelength,v0,n0,NULL);
 sx0=v0[0];
 sy0=v0[1];

 if ( out_center_x != NULL )	*out_center_x=sx0;
 if ( out_center_y != NULL )	*out_center_y=sy0;

 cnt=0;
 for ( i=1 ; i<=nrings ; i++ )
  {	nsec=i*6;
	for ( j=0 ; j<nsec ; j++ )
	 {	double	x,y,r,w;
		r=aperture_radius*(double)i/(double)nrings;
		w=2*M_PI*j/(double)nsec;
		x=r*cos(w);
		y=r*sin(w);
		v0[0]=x;
		v0[1]=y;
		v0[2]=zstart-0.0001;
		n0[0]=angle_nx;
		n0[1]=angle_ny;
		n0[2]=angle_nz;
		optcalc_ray_trace(opt,wavelength,v0,n0,NULL);
		if ( cnt<out_max )
		 {	out_xy[cnt*2+0]=(v0[0]-sx0)/pixel_scale;
			out_xy[cnt*2+1]=(v0[1]-sy0)/pixel_scale;
		 }
		cnt++;
	 }
  }

 if ( out_n != NULL )	*out_n=cnt;

 return(0);
}

/*****************************************************************************/
/* gropt_single_raytrace: trace a single ray through the system              */
/*****************************************************************************/

int gropt_single_raytrace(struct optics *opt,double wavelength,
	double x0,double y0,double z0,
	double nx,double ny,double nz,
	double *out_points,int out_max,int *out_n)
{
 vector		v0,n0;
 struct raytrace	rt;
 int		i;

 optcalc_glass_refraction_precompute(opt,wavelength);
 optcalc_raytrace_reset(&rt);

 v0[0]=x0; v0[1]=y0; v0[2]=z0;
 n0[0]=nx; n0[1]=ny; n0[2]=nz;

 optcalc_ray_trace(opt,wavelength,v0,n0,&rt);

 if ( out_n != NULL )	*out_n=rt.rt_npoint;
 for ( i=0 ; i<rt.rt_npoint && i<out_max ; i++ )
  {	out_points[i*3+0]=rt.rt_points[i][0];
	out_points[i*3+1]=rt.rt_points[i][1];
	out_points[i*3+2]=rt.rt_points[i][2];
  }

 optcalc_raytrace_free(&rt);

 return(0);
}

/*****************************************************************************/
/* gropt_export_scad: export OpenSCAD 3D model as string                     */
/*****************************************************************************/

int gropt_export_scad(struct optics *opt,
	char **out_str,size_t *out_len)
{
 int		i;
 gropt_strbuf	sb;
 static char *openscad_color[]={"red","green","blue","magenta","darkcyan","brown","grey"};

 gropt_strbuf_init(&sb);

 gropt_strbuf_appendf(&sb,"/* generated by pyfitsh gropt */\n");
 for ( i=0 ; i<opt->opt_nlens ; i++ )
  {	struct	lens	*l;
	struct	surface	*s1,*s2;
	int		j,jmax;
	double		z0,z1;
	l=&opt->opt_lenses[i];
	z0=l->l_offset-l->l_thickness/2.0;
	z1=l->l_offset+l->l_thickness/2.0;
	jmax=32;
	gropt_strbuf_appendf(&sb,"color(\"%s\")rotate_extrude($fn=200)polygon([[%g,%g]",openscad_color[i%7],0.0,z0);
	s1=&l->l_s1;
	s2=&l->l_s2;
	for ( j=1 ; j<=jmax ; j++ )
	 {	double	r,z;
		r=l->l_radius1*(double)j/(double)jmax;
		optcalc_surface_aspheric_eval(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r*r,&z);
		gropt_strbuf_appendf(&sb,",[%g,%g]",r,z0+z);
	 }
	if ( l->l_radius1<l->l_radius2 )
	 {	double	r,z;
		r=l->l_radius1;
		optcalc_surface_aspheric_eval(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r*r,&z);
		gropt_strbuf_appendf(&sb,",[%g,%g]",l->l_radius2,z0+z);
	 }
	else if ( l->l_radius2<l->l_radius1 )
	 {	double	r,z;
		r=l->l_radius2;
		optcalc_surface_aspheric_eval(s2->s_curvature,s2->s_conic,s2->s_nalpha,s2->s_alphas,r*r,&z);
		gropt_strbuf_appendf(&sb,",[%g,%g]",l->l_radius1,z1-z);
	 }
	for ( j=jmax ; 0<=j ; j-- )
	 {	double	r,z;
		r=l->l_radius2*(double)j/(double)jmax;
		optcalc_surface_aspheric_eval(s2->s_curvature,s2->s_conic,s2->s_nalpha,s2->s_alphas,r*r,&z);
		gropt_strbuf_appendf(&sb,",[%g,%g]",r,z1-z);
	 }
	gropt_strbuf_appendf(&sb,"]);\n");
  }

 *out_str=sb.data;
 if ( out_len != NULL )	*out_len=sb.len;

 return(0);
}

/*****************************************************************************/
/* gropt_export_eps: export EPS planar diagram as string                     */
/*****************************************************************************/

int gropt_export_eps(struct optics *opt,
	double aperture_radius,int nrings,
	double angle_nx,double angle_ny,double zstart,
	double wavelength,double pixel_scale,
	char **out_str,size_t *out_len)
{
 int		i;
 double		ox,oy,scale,rmax,zmin,zmax,sx,sy,bndx,bndy;
 double		angle_nz;
 gropt_strbuf	sb;

 angle_nz=sqrt(1.0-(angle_nx*angle_nx+angle_ny*angle_ny));

 gropt_strbuf_init(&sb);

 rmax=0.0;
 zmin=zmax=opt->opt_z_focal;
 for ( i=0 ; i<opt->opt_nlens ; i++ )
  {	struct	lens	*l;
	double	z0,z1;
	l=&opt->opt_lenses[i];
	if ( rmax<l->l_radius1 )	rmax=l->l_radius1;
	if ( rmax<l->l_radius2 )	rmax=l->l_radius2;
	z0=l->l_offset-l->l_thickness/2.0;
	z1=l->l_offset+l->l_thickness/2.0;
	if ( z0<zmin )	zmin=z0;
	if ( zmax<z1 )	zmax=z1;
  }

 zmin=floor(zmin);
 zmax=1+floor(zmax);
 rmax=1+floor(rmax);

 scale=2.0;
 bndx=10.0;
 bndy=10.0;
 ox=scale*(-zmin+bndx);
 sx=scale*(zmax-zmin+2*bndx);
 oy=scale*(rmax+bndy);
 sy=2*oy;

 gropt_strbuf_appendf(&sb,"%%!PS-Adobe-2.0\n");
 gropt_strbuf_appendf(&sb,"%%%%Title: pyfitsh gropt\n");
 gropt_strbuf_appendf(&sb,"%%%%Creator: pyfitsh gropt\n");
 gropt_strbuf_appendf(&sb,"%%%%BoundingBox: 0 0 %g %g\n",sx,sy);

 for ( i=0 ; i<opt->opt_nlens ; i++ )
  {	struct	lens	*l;
	struct	surface	*s0,*s1;
	double		z0,z1,r0,r1;
	double		c0,dz0,zh0,zf0;
	double		c1,dz1,zh1,zf1;
	l=&opt->opt_lenses[i];
	z0=l->l_offset-l->l_thickness/2.0;
	z1=l->l_offset+l->l_thickness/2.0;
	gropt_strbuf_appendf(&sb,"%g %g moveto\n",ox+scale*z0,oy);
	s0=&l->l_s1;
	s1=&l->l_s2;
	r0=l->l_radius1;
	r1=l->l_radius2;
	optcalc_surface_aspheric_eval(s0->s_curvature,s0->s_conic,s0->s_nalpha,s0->s_alphas,r0*r0/4,&zh0);
	optcalc_surface_aspheric_eval(s0->s_curvature,s0->s_conic,s0->s_nalpha,s0->s_alphas,r0*r0  ,&zf0);
	optcalc_surface_aspheric_diff(s0->s_curvature,s0->s_conic,s0->s_nalpha,s0->s_alphas,r0*r0  ,&dz0);
	dz0=2*r0*dz0;
	if ( r0*dz0 != 0 )	c0=8.0*(0.5*zf0-zh0)/(3.0*r0*dz0);
	else			c0=1.0/3.0;
	optcalc_surface_aspheric_eval(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r1*r1/4,&zh1);
	optcalc_surface_aspheric_eval(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r1*r1  ,&zf1);
	optcalc_surface_aspheric_diff(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r1*r1  ,&dz1);
	dz1=2*r1*dz1;
	if ( r1*dz1 != 0 )	c1=8.0*(0.5*zf1-zh1)/(3.0*r1*dz1);
	else			c1=1.0/3.0;
	gropt_strbuf_appendf(&sb,"%g %g %g %g %g %g curveto\n",
		ox+scale*z0,oy+scale*(c0*r0),
		ox+scale*(z0+zf0-c0*r0*dz0),oy+scale*(1-c0)*r0,
		ox+scale*(z0+zf0),oy+scale*r0);
	if ( r0<r1 )
		gropt_strbuf_appendf(&sb,"%g %g lineto %g %g lineto\n",
			ox+scale*(z0+zf0),oy+scale*r1,
			ox+scale*(z1-zf1),oy+scale*r1);
	else if ( r1<r0 )
		gropt_strbuf_appendf(&sb,"%g %g lineto %g %g lineto\n",
			ox+scale*(z1-zf1),oy+scale*r0,
			ox+scale*(z1-zf1),oy+scale*r1);
	else
		gropt_strbuf_appendf(&sb,"%g %g lineto\n",ox+scale*(z1-zf1),oy+scale*r0);
	gropt_strbuf_appendf(&sb,"%g %g %g %g %g %g curveto\n",
		ox+scale*(z1-zf1+c1*r1*dz1),oy+scale*(1-c1)*r1,
		ox+scale*z1,oy+scale*(c1*r1),
		ox+scale*z1,oy);
	gropt_strbuf_appendf(&sb,"%g %g %g %g %g %g curveto\n",
		ox+scale*z1,oy-scale*(c1*r1),
		ox+scale*(z1-zf1+c1*r1*dz1),oy-scale*(1-c1)*r1,
		ox+scale*(z1-zf1),oy-scale*r1);
	if ( r0<r1 )
		gropt_strbuf_appendf(&sb,"%g %g lineto %g %g lineto\n",
			ox+scale*(z0+zf0),oy-scale*r1,
			ox+scale*(z0+zf0),oy-scale*r0);
	else if ( r1<r0 )
		gropt_strbuf_appendf(&sb,"%g %g lineto %g %g lineto\n",
			ox+scale*(z1-zf1),oy-scale*r0,
			ox+scale*(z0+zf0),oy-scale*r0);
	else
		gropt_strbuf_appendf(&sb,"%g %g lineto\n",
			ox+scale*(z0+zf0),oy-scale*r0);
	gropt_strbuf_appendf(&sb,"%g %g %g %g %g %g curveto\n",
		ox+scale*(z0+zf0-c0*r0*dz0),oy-scale*(1-c0)*r0,
		ox+scale*z0,oy-scale*(c0*r0),
		ox+scale*z0,oy);
	gropt_strbuf_appendf(&sb,"gsave 0.9 setgray fill grestore\n");
	gropt_strbuf_appendf(&sb,"1 setlinewidth stroke\n");
  }

 if ( 0<nrings && 0.0<aperture_radius )
  {	optcalc_glass_refraction_precompute(opt,wavelength);
	for ( i=-nrings ; i<=nrings ; i++ )
	 {	double	r;
		vector	v0,n0;
		struct	raytrace	rt;
		int	j;
		if ( nrings )
			r=aperture_radius*(double)i/(double)nrings;
		else
			r=0.0;
		v0[0]=0.0;
		v0[1]=r;
		v0[2]=zstart-0.0001;
		n0[0]=angle_nx;
		n0[1]=angle_ny;
		n0[2]=angle_nz;
		optcalc_raytrace_reset(&rt);
		optcalc_ray_trace(opt,wavelength,v0,n0,&rt);
		for ( j=0 ; j<rt.rt_npoint ; j++ )
		 {	gropt_strbuf_appendf(&sb,"%g %g %s\n",
				ox+scale*rt.rt_points[j][2],oy+scale*rt.rt_points[j][1],
				j?"lineto":"moveto");
		 }
		gropt_strbuf_appendf(&sb,"0.5 setlinewidth stroke\n");
		optcalc_raytrace_free(&rt);
	 }
  }

 *out_str=sb.data;
 if ( out_len != NULL )	*out_len=sb.len;

 return(0);
}

/*****************************************************************************/
/* gropt_compute_psf: point-spread function via triangle adaptive subdivision*/
/*****************************************************************************/

typedef struct triangle_s triangle;

typedef struct
 {	double	x;
	double	y;
 } triangle_point;

struct triangle_s
 {	triangle	*prev,*next;
	triangle_point	t_vertex[3];
	int		t_flag;
	double		t_area,t_proj;
 };

int gropt_compute_psf(struct optics *opt,
	double wavelength,double aperture_radius,
	double angle_nx,double angle_ny,double zstart,double pixel_scale,
	int psf_hsize,double *out_psf)
{
 triangle	*triangles,*t;
 int		i,iteration_depth,max_iteration;
 double		*polygon_buf,sx0,sy0;
 int		px,py,psf_size;
 vector		v0,n0;
 double		angle_nz;

 angle_nz=sqrt(1.0-(angle_nx*angle_nx+angle_ny*angle_ny));

 optcalc_glass_refraction_precompute(opt,wavelength);

 v0[0]=0.0;
 v0[1]=0.0;
 v0[2]=zstart-0.0001;
 n0[0]=angle_nx;
 n0[1]=angle_ny;
 n0[2]=angle_nz;
 optcalc_ray_trace(opt,wavelength,v0,n0,NULL);
 sx0=v0[0];
 sy0=v0[1];

 triangles=NULL;
 for ( i=0 ; i<6 ; i++ )
  {	double	w1,w2;
	t=list_new(triangle);
	t->t_vertex[0].x=0.0;
	t->t_vertex[0].y=0.0;
	w1=(double)(i+0)*M_PI/3.0;
	w2=(double)(i+1)*M_PI/3.0;
	t->t_vertex[1].x=aperture_radius*cos(w1);
	t->t_vertex[1].y=aperture_radius*sin(w1);
	t->t_vertex[2].x=aperture_radius*cos(w2);
	t->t_vertex[2].y=aperture_radius*sin(w2);
	t->t_flag=1;
	list_insert_first(triangles,t);
  }

 max_iteration=7;
 for ( iteration_depth=0 ; iteration_depth<max_iteration ; iteration_depth++ )
  {	for ( t=triangles ; t != NULL ; )
	 {	int	ii,cnt;
		if ( ! t->t_flag )
		 {	t=t->next;
			continue;
		 }
		cnt=0;
		for ( ii=0 ; ii<3 ; ii++ )
		 {	v0[0]=t->t_vertex[ii].x;
			v0[1]=t->t_vertex[ii].y;
			v0[2]=zstart-0.0001;
			n0[0]=angle_nx;
			n0[1]=angle_ny;
			n0[2]=angle_nz;
			if ( ! optcalc_ray_trace(opt,wavelength,v0,n0,NULL) )
				cnt++;
		 }
		if ( cnt==3 )
		 {	t->t_flag=0;
			t=t->next;
			continue;
		 }
		else if ( cnt==0 )
		 {	triangle	*n;
			n=t->next;
			list_remove(triangles,t);
			free(t);
			t=n;
		 }
		else
		 {	triangle	*n,*t1,*t2,*t3,*t4;
			triangle_point	v01,v12,v20;
			n=t->next;
			t1=list_new(triangle);
			t2=list_new(triangle);
			t3=list_new(triangle);
			t4=list_new(triangle);
			v01.x=(t->t_vertex[0].x+t->t_vertex[1].x)/2.0;
			v01.y=(t->t_vertex[0].y+t->t_vertex[1].y)/2.0;
			v12.x=(t->t_vertex[1].x+t->t_vertex[2].x)/2.0;
			v12.y=(t->t_vertex[1].y+t->t_vertex[2].y)/2.0;
			v20.x=(t->t_vertex[2].x+t->t_vertex[0].x)/2.0;
			v20.y=(t->t_vertex[2].y+t->t_vertex[0].y)/2.0;
			t1->t_vertex[0]=t->t_vertex[0];
			t1->t_vertex[1]=v01;
			t1->t_vertex[2]=v20;
			t1->t_flag=1;
			t2->t_vertex[0]=t->t_vertex[1];
			t2->t_vertex[1]=v12;
			t2->t_vertex[2]=v01;
			t2->t_flag=1;
			t3->t_vertex[0]=t->t_vertex[2];
			t3->t_vertex[1]=v20;
			t3->t_vertex[2]=v12;
			t3->t_flag=1;
			t4->t_vertex[0]=v01;
			t4->t_vertex[1]=v12;
			t4->t_vertex[2]=v20;
			t4->t_flag=1;
			list_remove(triangles,t);
			free(t);
			list_insert_first(triangles,t1);
			list_insert_first(triangles,t2);
			list_insert_first(triangles,t3);
			list_insert_first(triangles,t4);
			t=n;
		 }
	 }
  }

 for ( t=triangles ; t != NULL ; t=t->next )
  {	double	dx1,dy1,dx2,dy2;
	int	ii;
	if ( t->t_flag )	continue;
	dx1=t->t_vertex[1].x-t->t_vertex[0].x;
	dy1=t->t_vertex[1].y-t->t_vertex[0].y;
	dx2=t->t_vertex[2].x-t->t_vertex[0].x;
	dy2=t->t_vertex[2].y-t->t_vertex[0].y;
	t->t_area=fabs(dx1*dy2-dx2*dy1)/2.0;
	for ( ii=0 ; ii<3 ; ii++ )
	 {	v0[0]=t->t_vertex[ii].x;
		v0[1]=t->t_vertex[ii].y;
		v0[2]=zstart-0.0001;
		n0[0]=angle_nx;
		n0[1]=angle_ny;
		n0[2]=angle_nz;
		optcalc_ray_trace(opt,wavelength,v0,n0,NULL);
		t->t_vertex[ii].x=(v0[0]-sx0)/pixel_scale;
		t->t_vertex[ii].y=(v0[1]-sy0)/pixel_scale;
	 }
	dx1=t->t_vertex[1].x-t->t_vertex[0].x;
	dy1=t->t_vertex[1].y-t->t_vertex[0].y;
	dx2=t->t_vertex[2].x-t->t_vertex[0].x;
	dy2=t->t_vertex[2].y-t->t_vertex[0].y;
	t->t_proj=fabs(dx1*dy2-dx2*dy1)/2.0;
  }

 psf_size=2*psf_hsize+1;
 polygon_buf=(double *)malloc(sizeof(double)*2*16);

 for ( py=-psf_hsize ; py<=psf_hsize ; py++ )
  {	for ( px=-psf_hsize ; px<=psf_hsize ; px++ )
	 {	double	a;
		a=0.0;
		for ( t=triangles ; t != NULL ; t=t->next )
		 {	int	j,n;
			double	w;
			if ( t->t_flag )	continue;
			for ( j=0 ; j<3 ; j++ )
			 {	polygon_buf[2*j+0]=t->t_vertex[j].x;
				polygon_buf[2*j+1]=t->t_vertex[j].y;
			 }
			n=polygon_intersection_square(polygon_buf,3,(double)px-0.5,(double)py-0.5,1.0,1.0);
			if ( n<=0 )	continue;
			w=polygon_area(polygon_buf,n);
			a+=t->t_area*w/t->t_proj;
		 }
		out_psf[(psf_hsize+py)*psf_size+(psf_hsize+px)]=a;
	 }
  }

 free(polygon_buf);

 while ( triangles != NULL )
  {	t=triangles;
	list_remove(triangles,t);
	free(t);
  }

 optcalc_glass_refraction_precompute(opt,0.0);

 return(0);
}
