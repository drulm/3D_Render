
/*----------------------------------------------------------------
Phong, Gourand, Flat shading program (wireframe also)
for Linux X windows.

Coded in 3 days. Prototype & proof of concept only,
Not for production. 
----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#include <X11/IntrinsicP.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>
#include <X11/Core.h>

/*--------------------------------------------------------------------*/

int light();
double light2();
double light3();
double dist();
void plane();
void render_triangles();
void draw_lines_back();
void load_object();
void translate();
void get_edge2();
void rot_obj();
void plane_calc();
void device_trans_all();
int getcode();
void draw_lines();
void clip();
void phong();
void get_pnt();
void put_pnt();
void copy_all();
void gouraud();
void proj_trans_all();
void proj_point();
void build_index();
void close_all();
void get_point();
void get_edge();
void make_elipse();
void view_trans();
void view_trans_all();
void add_to_list();
void set_line();
void init();
double mag();
void set_point();
double to_rad();
double pw();
void cross();
void find_nuv();
void rotate();
void get_line();
void load_scene();
/*--------------------------------------------------------------------*/

#define MAX_COLORS	1024
#define LIGHT_MAX	10
#define NUM_OBJS	20
#define SCALEFAC	500.0
#define PI		3.14159
#define TESTING_X	0
#define	TESTING_Y	0
#define	TESTING_WIDTH	400
/*#define	TESTING_WIDTH	800 */
#define CLIP		100000.0
#define	TESTING_HEIGHT  TESTING_WIDTH
#define M		TESTING_WIDTH/2
#define	TESTING_BDWIDTH 1
#define MAX_ELEM	100000
/*--------------------------------------------------------------------*/

static char TESTING_DEFAULTFONT[] = "9x15";
static char TESTING_TITLE[] = "Superellipsoid";


struct line
{
	double x1, y1, z1, x2, y2, z2;		/* obj modeling only! */
};


struct point
{
	double x, y, z;					/* a point */
};


struct edge
{
	int p1, p2;						/* pointer into edge table */
};


struct poly
{
	int st, end;						/* st-end in edge table */
	double a, b, c, d;						/* plane equation */
};


struct objrange
{
	int st, end;				/* st-end of obj edges */
	double x, y, z;				/* location of obj */
	double rx, ry, rz;			/* rotation of obj */
	double ka;				/* ambient reflect */
	double kd;				/* diffuse reflect */
	double ks;				/* spect reflect */
	double kn;				/* shinyness */
	double dx, dy, dz;			/* move */
	double tx, ty, tz;			/* texture */
	double ttx, tty, ttz;			/* texture */
	double drx, dry, drz;			/* spin! */
};


struct light
{
	double x, y, z;				/* loc of light */
	double I;				/* intensity */
};


struct index
{
	int org;
	int poly;
};


struct point_norm
{
	double a, b, c;			/* plane equation */
	double l;
	int lp;					/* A POLYGON */
	int new;				/* points to rotated new pt */
};


/*--------------------------------------------------------------------*/
Display		*display;
Window		topwin;
Pixmap		toppix;
XFontStruct	*fontSt;
GC			topwinGc;
XEvent		event;

int			default_depth;
Visual		*default_visual;
static		char *name[] = { "Red", "Yellow", "Green" };
XColor		exact_defs[MAX_COLORS];
XColor		color;
Colormap	default_cmap;
int			ncolors = MAX_COLORS;
int			plane_masks[1];
int			colors[MAX_COLORS];
XVisualInfo visual_info;
int class;

XtAppContext app;
unsigned long NumColors;

double zbuff[TESTING_HEIGHT + 1][TESTING_WIDTH + 1];

struct line a[MAX_ELEM], b[MAX_ELEM], bo[MAX_ELEM], rot;
struct point VRP, N, V, X, n, u, v;
struct point PRP, Wmin, Wmax, Cmin, Cmax;
double XSCALE, YSCALE;
int seg, sv;

double AMB;
struct light lg[LIGHT_MAX];
struct point pt[MAX_ELEM], ptn[MAX_ELEM];
struct point_norm nm[MAX_ELEM];
struct index idx[MAX_ELEM];
short int roted[MAX_ELEM];
struct edge et[MAX_ELEM], etn[MAX_ELEM];
struct poly poly[MAX_ELEM / 3];
struct objrange obj[NUM_OBJS];	/* pointer into edge table to list of objs */
int pntcnt, pntcnt2, edgecnt, edgecnt2, objcnt, plycnt, lgcnt;

int	screen_num;
unsigned long	fgColor, bgColor, bdColor;

int i;
int DIV;
double x, y, z;
double xp, yp, zp;
double rx, ry, rz;
double s1, s2;
double g, t;
double ax, ay, az;
double dx, dy, dz;


/*-----------------------------------------------------------------------*/
/*--------------------------------------- MAIN --------------------------*/
/*-----------------------------------------------------------------------*/
main(argc, argv)
int argc;
char *argv[];
{
	/*	Initialization - connect to server and get screen	*/
	init();
	pntcnt = -1;
	edgecnt = -1;
	objcnt = -1;
	plycnt = -1;
	lgcnt = -1;

	set_point(&VRP, 0.0, 0.0, 0.0);
	set_point(&N, 0.0, 0.0, 1.0);
	set_point(&V, 0.0, 1.0, 0.0);

	find_nuv(&VRP, &N, &V, &n, &u, &v);

	set_point(&PRP, 0.0, 0.0, 1000.0);
	set_point(&Wmin, -1000.0, -1000.0, 0.0);
	set_point(&Wmax, 1000.0, 1000.0, 0.0);
	set_point(&Cmin, -CLIP, -CLIP, -CLIP);
	set_point(&Cmax, CLIP, CLIP, 990.0);

	XSCALE = ((double)TESTING_WIDTH - TESTING_X) / (Wmax.x - Wmin.x);
	YSCALE = ((double)TESTING_HEIGHT - TESTING_Y) / (Wmax.y - Wmin.y);

	load_scene();
	//load_object();
	//make_elipse();

	/*--------------------- DO the object rotation animation ----------------*/
	printf("** ANIMATE **\n\n");
	ax = 1.0; ay = 1.0; az = 1.0;
	dx = 0.0;
	dy = 0.0;
	dz = 0.0;

	copy_all();
	build_index();
	printf("SIZES (%d,%d)\n", pntcnt, pntcnt2);
	XNextEvent(display, &event);

	while (1)
	{
		copy_all();			/* FIX!!! TO MAKE JUST LINE PRIM. */
		for (i = 0; i <= objcnt; i++)
		{
			rot_obj(i);
			obj[i].x += obj[i].dx;
			obj[i].y += obj[i].dy;
			obj[i].z += obj[i].dz;
			obj[i].rx += obj[i].drx;
			obj[i].ry += obj[i].dry;
			obj[i].rz += obj[i].drz;
		}

		view_trans_all(&VRP, &n, &u, &v);
		proj_trans_all(&PRP, &Wmin, &Wmax);

		clip();
		device_trans_all();
		plane_calc();

		XSetForeground(display, topwinGc, bgColor);
		XFillRectangle(display, toppix, topwinGc, 0, 0,
			TESTING_WIDTH, TESTING_HEIGHT);
		XSetForeground(display, topwinGc, fgColor);

		draw_lines_back();

		/*		draw_lines_back();

		XCopyArea(display,toppix,topwin,topwinGc,
		0,0, TESTING_WIDTH, TESTING_HEIGHT, 0, 0);

		render_triangles();

		XCopyArea(display,toppix,topwin,topwinGc,
		0,0, TESTING_WIDTH, TESTING_HEIGHT, 0, 0);

		gouraud();

		XCopyArea(display,toppix,topwin,topwinGc,
		0,0, TESTING_WIDTH, TESTING_HEIGHT, 0, 0);

		*/
		//phong();

		XCopyArea(display, toppix, topwin, topwinGc,
			0, 0, TESTING_WIDTH, TESTING_HEIGHT, 0, 0);

		XNextEvent(display, &event);
	}
	close_all();
}


/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/// <summary>
/// Build_indexes this instance.
/// </summary>
void build_index()
{
	int st, end, i, j, p1, p2;
	int pc2;

	fprintf(stderr, "	*PRECOMPUTATION*\n");

	for (i = 0; i <= plycnt; i++)
	{				/* find polygon for each pt */
		if (i % 100 == 0)
		{
			fprintf(stderr, "precomputed %d/%d polygons\n", i, plycnt);
		}

		st = poly[i].st;
		end = poly[i].end;

		for (j = st; j <= end; j++)
		{
			p1 = etn[j].p1; 
			p2 = etn[j].p2;

			idx[p1].poly = i; 
			idx[p2].poly = i;
		}
	}

	pc2 = -1;

	/* loop over lines and pull out NEW points...set new # of pts */
	for (i = 0; i <= edgecnt; i++)
	{
		p1 = et[i].p1; 
		p2 = et[i].p2;

		pc2++;
		idx[pc2].org = p1;
		nm[p1].new = pc2;

		pc2++;
		idx[pc2].org = p2;
		nm[p2].new = pc2;
	}
}

/*--------------------------------------------------------------------*/
double light3(p, nx, ny, nz, x, y, z)
int p;			/* the triangle lit */
double nx, ny, nz, x, y, z;
{
	int i, st, end, pt, pp[20], pps, j, in, ob;
	int new;
	struct point vec, vc2;
	double H, vx, vy, vz, lx, ly, lz, hx, hy, hz, LIGHT, LT, m1;

	st = poly[p].st;
	end = poly[p].end;

	vx = vy = vz = 0;

	pps = -1;
	LIGHT = 0.0;

	for (i = 0; i <= objcnt; i++)
	{
		if (obj[i].st <= p && obj[i].end >= p)
		{
			ob = i;
			break;
		}
	}

	if (obj[ob].tx > 0.0)
	{
		nx += sin(x / obj[ob].tx) / obj[ob].ttx;
		ny += sin(y / obj[ob].ty) / obj[ob].tty;
		nz += sin(z / obj[ob].tz) / obj[ob].ttz;
	}

	LIGHT += obj[ob].ka * AMB;
	for (i = 0; i <= lgcnt; i++)
	{
		lx = lg[i].x - x; 
		ly = lg[i].y - y; 
		lz = lg[i].z - z;

		vec.x = lx; 
		vec.y = ly; 
		vec.z = lz;

		m1 = mag(&vec);
		lx /= m1; 
		ly /= m1; 
		lz /= m1;

		vec.x = nx; 
		vec.y = ny; 
		vec.z = nz;

		m1 = mag(&vec);
		vec.x /= m1; vec.y /= m1; vec.z /= m1;
		LT = obj[ob].kd * lg[i].I * (lx*vec.x + ly*vec.y + lz*vec.z);

		if (LT > 0)
		{
			LIGHT += LT;
		}

		vx = PRP.x - x; vy = PRP.y - y; vz = PRP.z - z;
		vc2.x = vx; vc2.y = vy; vc2.z = vz;

		m1 = mag(&vc2);
		vx /= m1; vy /= m1; vz /= m1;
		hx = lx + vx; hy = ly + vy; hz = lz + vz;
		vc2.x = hx; vc2.y = hy; vc2.z = hz;

		m1 = mag(&vc2);
		hx /= m1; hy /= m1; hz /= m1;
		H = hx*vec.x + hy*vec.y + hz*vec.z;

		if (H > 0.0)
		{
			LT = obj[ob].ks * lg[i].I * pow(H, obj[ob].kn);
		}
		else
		{
			LT = -1.0;
		}

		if (LT > 0)
		{
			LIGHT += LT;
		}
	}

	if (LIGHT > 1.0)
	{
		LIGHT = 1.0;
	}

	return((LIGHT*((double)NumColors - 1.0)));
}

/*--------------------------------------------------------------------*/
double light2(p, v)
int p, v;			/* the triangle lit */
{
	int i, st, end, pt, pp[20], pps, j, in, ob;
	int new;
	struct point vec, vc2;
	double H, vx, vy, vz, lx, ly, lz, hx, hy, hz, LIGHT, LT, m1;

	st = poly[p].st;
	end = poly[p].end;

	vx = vy = vz = 0;
	pps = -1;
	LIGHT = 0.0;

	for (i = 0; i <= objcnt; i++)
	{
		if (obj[i].st <= p && obj[i].end >= p)
		{
			ob = i;
			break;
		}
	}

	LIGHT += obj[ob].ka * AMB;

	for (i = 0; i <= lgcnt; i++)
	{
		new = nm[v].new;
		lx = lg[i].x - ptn[new].x;
		ly = lg[i].y - ptn[new].y;
		lz = lg[i].z - ptn[new].z;

		vec.x = lx;
		vec.y = ly;
		vec.z = lz;

		m1 = mag(&vec);
		lx /= m1;
		ly /= m1;
		lz /= m1;

		vec.x = nm[v].a;
		vec.y = nm[v].b;
		vec.z = nm[v].c;

		m1 = mag(&vec);
		vec.x /= m1;
		vec.y /= m1;
		vec.z /= m1;

		LT = obj[ob].kd * lg[i].I * (lx*vec.x + ly*vec.y + lz*vec.z);

		if (LT > 0)
		{
			LIGHT += LT;
		}

		vx = 0.0 - ptn[new].x;
		vy = 0.0 - ptn[new].y;
		vz = 1000.0 - ptn[new].z;
		vc2.x = vx;
		vc2.y = vy;
		vc2.z = vz;

		m1 = mag(&vc2);
		vx /= m1;
		vy /= m1;
		vz /= m1;

		hx = lx + vx;
		hy = ly + vy;
		hz = lz + vz;

		vc2.x = hx;
		vc2.y = hy;
		vc2.z = hz;

		m1 = mag(&vc2);
		hx /= m1;
		hy /= m1;
		hz /= m1;

		H = hx*vec.x + hy*vec.y + hz*vec.z;
		if (H > 0.0)
		{
			LT = obj[ob].ks * lg[i].I * pow(H, obj[ob].kn);
		}
		else
		{
			LT = -1.0;
		}

		if (LT > 0)
		{
			LIGHT += LT;
		}
	}

	if (LIGHT > 1.0)
	{
		LIGHT = 1.0;
	}

	return((LIGHT*((double)NumColors - 1.0)));
}


/*--------------------------------------------------------------------*/
int light(p)
int p;			/* the triangle lit */
{
	int i, st, end, pt, pp[20], pps, j, in, ob;
	struct point vec;
	double vx, vy, vz, lx, ly, lz, LIGHT, LT, m1;

	st = poly[p].st;
	end = poly[p].end;

	vx = vy = vz = 0;

	pps = -1;
	LIGHT = 0.0;

	for (i = 0; i <= objcnt; i++)
		if (obj[i].st <= p && obj[i].end >= p)
		{
			ob = i;
			break;
		}

	LIGHT += obj[ob].ka * AMB;

	for (i = st; i <= end; i++)
	{
		pt = etn[i].p1;
		in = 0;

		for (j = 0; j <= pps; j++)
		{
			if (pt == pp[j])
			{
				in = 1;
			}
		}

		if (!in)
		{
			pps++;
			pp[pps] = pt;

			vx += ptn[pt].x;
			vy += ptn[pt].y;
			vz += ptn[pt].z;
		}

		pt = etn[i].p2;
		in = 0;

		for (j = 0; j <= pps; j++)
		{
			if (pt == pp[j])
			{
				in = 1;
			}
		}

		if (!in)
		{
			pps++;
			pp[pps] = pt;

			vx += ptn[pt].x;  
			vy += ptn[pt].y;  
			vz += ptn[pt].z;
		}
	}

	vx /= (pps + 1); 
	vy /= (pps + 1); 
	vz /= (pps + 1);

	for (i = 0; i <= lgcnt; i++)
	{
		lx = lg[i].x - vx; 
		ly = lg[i].y - vy; 
		lz = lg[i].z - vz;

		vec.x = lx; 
		vec.y = ly; 
		vec.z = lz;

		m1 = mag(&vec);

		lx /= m1; 
		ly /= m1; 
		lz /= m1;

		vec.x = poly[p].a; 
		vec.y = poly[p].b; 
		vec.z = poly[p].c;

		m1 = mag(&vec);
		LT = obj[ob].kd * lg[i].I * (lx*poly[p].a / m1 + ly*poly[p].b / m1 + lz*poly[p].c / m1);

		if (LT > 0)
		{
			LIGHT += LT;
		}
	}

	if (LIGHT > 1.0)
	{
		LIGHT = 1.0;
	}

	return(exact_defs[(int)(LIGHT*(NumColors - 1))].pixel);
}


/*--------------------------------------------------------------------*/
void load_scene()
{
	FILE *fs;
	char name[80], com;

	printf("enter a scene file >");
	gets(name);
	fs = fopen(name, "r");

	if (fs == NULL)
	{
		exit(-1);
	}

	com = '*';
	while (!feof(fs) && com != 'Q')
	{
		fscanf(fs, "%c\n", &com);

		if (com == 'A')
		{
			fscanf(fs, "%lf\n", &AMB);
			printf("ambient light =%lf\n", AMB);
		}
		else if (com == 'L')
		{
			lgcnt++;
			fscanf(fs, "%lf%lf%lf%lf\n",
				&lg[lgcnt].x, &lg[lgcnt].y, &lg[lgcnt].z, &lg[lgcnt].I);
			printf("light %lf\n", lg[lgcnt].I);
		}
		else if (com == 'O')
		{
			fscanf(fs, "%s\n", name);
			printf("object: %s\n", name);
			load_object(name);
			fscanf(fs, "%lf%lf%lf\n"
				, &obj[objcnt].x, &obj[objcnt].y, &obj[objcnt].z);
			fscanf(fs, "%lf%lf%lf\n"
				, &obj[objcnt].rx, &obj[objcnt].ry, &obj[objcnt].rz);
			fscanf(fs, "%lf%lf%lf%lf\n"
				, &obj[objcnt].ka, &obj[objcnt].kd, &obj[objcnt].ks,
				&obj[objcnt].kn);
			fscanf(fs, "%lf%lf%lf%lf%lf%lf\n"
				, &obj[objcnt].tx, &obj[objcnt].ty, &obj[objcnt].tz
				, &obj[objcnt].ttx, &obj[objcnt].tty, &obj[objcnt].ttz);
			fscanf(fs, "%lf%lf%lf\n"
				, &obj[objcnt].dx, &obj[objcnt].dy, &obj[objcnt].dz);
			fscanf(fs, "%lf%lf%lf\n"
				, &obj[objcnt].drx, &obj[objcnt].dry, &obj[objcnt].drz);
		}
		else
		{
			printf("ERROR in SCENE FILE!\n");
			exit(-2);
		}
	}
}


/*--------------------------------------------------------------------*/
void phong()
{
	double x1[5], y1[5], z1[5];
	double x2[5], y2[5], z2[5];
	int e[5], ll, rl;
	double minx, maxx, miny, maxy, xp1, xp2, m, b, xx1, xx2, yy1, yy2;
	double m1, dep, IL, IR, IP, C1, C2;
	struct point vec, LN, RN, RP, llp, rrp, ppp;
	int i, j, k, st, end, c, lp, p1, p2;

	fprintf(stderr, "----PHONG SHADE----\n");
	for (k = 0; k <= pntcnt; k++)
	{
		nm[k].a = 0.0;
		nm[k].b = 0.0;
		nm[k].c = 0.0;
	}
	fprintf(stderr, "	*computing normals*\n");
	for (k = 0; k <= pntcnt2; k++)
	{
		i = idx[k].poly;
		j = idx[k].org;

		nm[j].a += poly[i].a;
		nm[j].b += poly[i].b;
		nm[j].c += poly[i].c;

		nm[j].lp = i;
	}
	fprintf(stderr, "	*lighting points*\n");
	for (k = 0; k <= pntcnt; k++)
	{
		vec.x = nm[k].a; vec.y = nm[k].b; vec.z = nm[k].c;

		m1 = mag(&vec);
		nm[k].a /= m1; nm[k].b /= m1; nm[k].c /= m1;
	}

	for (j = 0; j <= TESTING_HEIGHT; j++)
		for (k = 0; k <= TESTING_WIDTH; k++) zbuff[j][k] = -9999999999.0;

	for (i = 0; i <= plycnt; i++)
	{
		if (i % 100 == 0) fprintf(stderr, "rendered %d/%d polygons\n", i, plycnt);

		st = poly[i].st; 
		end = poly[i].end;

		c = 0;
		miny = (double)TESTING_WIDTH;
		maxy = 0.0;
		for (j = st; j <= end; j++)
			if (etn[j].p1 != -1)
			{
				e[c] = j;
				get_edge2(j, &x1[c], &y1[c], &z1[c], &x2[c], &y2[c], &z2[c]);

				if (y1[c]<miny) miny = y1[c];
				if (y2[c]<miny) miny = y2[c];
				if (y1[c]>maxy) maxy = y1[c];
				if (y2[c]>maxy) maxy = y2[c];

				c++;
			}

		for (k = floor(miny) + 1; k <= ceil(maxy) - 1; k++)
		{
			minx = (double)TESTING_WIDTH;
			maxx = 0.0;

			ll = -1;
			rl = -1;

			for (j = 0; j<c; j++)
			{
				if (x1[j]<x2[j]) { xx1 = x1[j]; xx2 = x2[j]; }
				else { xx1 = x2[j]; xx2 = x1[j]; }

				if (y1[j]<y2[j]) { yy1 = y1[j]; yy2 = y2[j]; }
				else { yy1 = y2[j]; yy2 = y1[j]; }

				if (x2[j] != x1[j]) m = (y2[j] - y1[j]) / (x2[j] - x1[j]);
				else m = 0.0000000000000001;

				b = y1[j] - m*x1[j];

				xp1 = ((double)k - b) / m;

				if (xp1<minx && xp1>xx1 &&
					(double)k>yy1 && (double)k<yy2) 
				{
					minx = xp1;  
					ll = e[j];
				}
				if (xp1>maxx && xp1<xx2 &&
					(double)k>yy1 && (double)k<yy2)
				{
					maxx = xp1;  
					rl = e[j];
				}
			}

			if (c == 3 && (int)minx>0 && (int)maxx<TESTING_WIDTH && ll != -1 && rl != -1)
			{
				p1 = etn[ll].p1;  
				p2 = etn[ll].p2;

				if (ptn[p1].y != ptn[p2].y)
				{
					C1 = ((double)k - ptn[p2].y) / (ptn[p1].y - ptn[p2].y);
					C2 = (ptn[p1].y - k) / (ptn[p1].y - ptn[p2].y);

					LN.x = C1*nm[idx[p1].org].a + C2*nm[idx[p2].org].a;
					LN.y = C1*nm[idx[p1].org].b + C2*nm[idx[p2].org].b;
					LN.z = C1*nm[idx[p1].org].c + C2*nm[idx[p2].org].c;

					llp.x = C1*ptn[p1].x + C2*ptn[p2].x;
					llp.y = C1*ptn[p1].y + C2*ptn[p2].y;
					llp.z = C1*ptn[p1].z + C2*ptn[p2].z;
				}
			    else
				{
					LN.x = nm[idx[p1].org].a;
					LN.y = nm[idx[p1].org].b;
					LN.z = nm[idx[p1].org].c;

					llp.x = ptn[p1].x;
					llp.y = ptn[p1].y;
					llp.z = ptn[p1].z;
				}

			    p1 = etn[rl].p1;  
				p2 = etn[rl].p2;

				if (ptn[p1].y != ptn[p2].y)
				{
					C1 = ((double)k - ptn[p2].y) / (ptn[p1].y - ptn[p2].y);
					C2 = (ptn[p1].y - k) / (ptn[p1].y - ptn[p2].y);

					RN.x = C1*nm[idx[p1].org].a + C2*nm[idx[p2].org].a;
					RN.y = C1*nm[idx[p1].org].b + C2*nm[idx[p2].org].b;
					RN.z = C1*nm[idx[p1].org].c + C2*nm[idx[p2].org].c;

					rrp.x = C1*ptn[p1].x + C2*ptn[p2].x;
					rrp.y = C1*ptn[p1].y + C2*ptn[p2].y;
					rrp.z = C1*ptn[p1].z + C2*ptn[p2].z;
				}
			    else
				{
					RN.x = nm[idx[p1].org].a;
					RN.y = nm[idx[p1].org].b;
					RN.z = nm[idx[p1].org].c;

					rrp.x = ptn[p1].x;
					rrp.y = ptn[p1].y;
					rrp.z = ptn[p1].z;
				}

			    for (j = floor(minx) + 1; j <= ceil(maxx) - 1; j++)
				{
					dep = (-poly[i].a*(double)j - poly[i].b*(double)k - poly[i].d) / poly[i].c;
					
					if (dep>zbuff[j][k])
					{
						if (maxx != minx)
						{
							C1 = (maxx - (double)j) / (maxx - minx);
							C2 = ((double)j - minx) / (maxx - minx);

							RP.x = C1*LN.x + C2*RN.x;
							RP.y = C1*LN.y + C2*RN.y;
							RP.z = C1*LN.z + C2*RN.z;

							ppp.x = C1*llp.x + C2*rrp.x;
							ppp.y = C1*llp.y + C2*rrp.y;
							ppp.z = C1*llp.z + C2*rrp.z;
						}
						else
						{
							RP.x = RN.x;
							RP.y = RN.y;
							RP.z = RN.z;

							ppp.x = rrp.x;
							ppp.y = rrp.y;
							ppp.z = rrp.z;
						}

						IP = light3(i, RP.x, RP.y, RP.z, ppp.x, ppp.y, ppp.z);

						if (IP>255.0) IP = 255.0; else if (IP<0.0) IP = 0.0;

						XSetForeground(display, topwinGc, exact_defs[(int)IP].pixel);

						XDrawPoint(display, toppix, topwinGc, j, k);

						zbuff[j][k] = dep;
					}
				}
			}
		}
	}
}

/*--------------------------------------------------------------------*/
void gouraud()
{
	double x1[5], y1[5], z1[5];
	double x2[5], y2[5], z2[5];
	int e[5], ll, rl;
	double minx, maxx, miny, maxy, xp1, xp2, m, b, xx1, xx2, yy1, yy2;
	double m1, dep, IL, IR, IP;
	struct point vec;
	int i, j, k, st, end, c, lp, p1, p2;

	fprintf(stderr, "----GOURAUD SHADE----\n");
	for (k = 0; k <= pntcnt; k++)
	{
		nm[k].a = 0.0;
		nm[k].b = 0.0; 
		nm[k].c = 0.0;
	}
	fprintf(stderr, "	*computing normals*\n");
	for (k = 0; k <= pntcnt2; k++)
	{
		i = idx[k].poly;
		j = idx[k].org;

		nm[j].a += poly[i].a;
		nm[j].b += poly[i].b;
		nm[j].c += poly[i].c;

		nm[j].lp = i;
	}
	fprintf(stderr, "	*lighting points*\n");
	for (k = 0; k <= pntcnt; k++)
	{
		vec.x = nm[k].a; 
		vec.y = nm[k].b; 
		vec.z = nm[k].c;

		m1 = mag(&vec);
		nm[k].a /= m1; 
		nm[k].b /= m1; 
		nm[k].c /= m1;

		nm[k].l = light2(nm[k].lp, k);
	}

	for (j = 0; j <= TESTING_HEIGHT; j++)
		for (k = 0; k <= TESTING_WIDTH; k++) 
			zbuff[j][k] = -9999999999.0;

	for (i = 0; i <= plycnt; i++)
	{
		if (i % 100 == 0) fprintf(stderr, "rendered %d/%d polygons\n", i, plycnt);

		st = poly[i].st; 
		end = poly[i].end;

		c = 0;
		miny = (double)TESTING_WIDTH;
		maxy = 0.0;

		for (j = st; j <= end; j++)
			if (etn[j].p1 != -1)
			{
				e[c] = j;
				get_edge2(j, &x1[c], &y1[c], &z1[c], &x2[c], &y2[c], &z2[c]);
				
				if (y1[c]<miny) miny = y1[c];
				if (y2[c]<miny) miny = y2[c];
				if (y1[c]>maxy) maxy = y1[c];
				if (y2[c]>maxy) maxy = y2[c];

				c++;
			}
		for (k = floor(miny) + 1; k <= ceil(maxy) - 1; k++)
		{
			minx = (double)TESTING_WIDTH;
			maxx = 0.0;

			ll = -1; 
			rl = -1;

			for (j = 0; j<c; j++)
			{
				if (x1[j]<x2[j])
					{ xx1 = x1[j]; xx2 = x2[j]; }
				else 
					{ xx1 = x2[j]; xx2 = x1[j]; }

				if (y1[j]<y2[j]) 
					{ yy1 = y1[j]; yy2 = y2[j]; }
				else 
					{ yy1 = y2[j]; yy2 = y1[j]; }
				
				if (x2[j] != x1[j]) 
					m = (y2[j] - y1[j]) / (x2[j] - x1[j]);
				else 
					m = 0.0000000000000001;

				b = y1[j] - m*x1[j];
				xp1 = ((double)k - b) / m;

				if (xp1<minx && xp1>xx1 && (double)k>yy1 && (double)k<yy2) 
				{
					minx = xp1;  ll = e[j];
				}

				if (xp1>maxx && xp1<xx2 && (double)k>yy1 && (double)k<yy2) 
				{
					maxx = xp1;  rl = e[j];
				}
			}

			if (c == 3 && (int)minx>0 && (int)maxx<TESTING_WIDTH
				&& ll != -1 && rl != -1)
			{
				p1 = etn[ll].p1;  
				p2 = etn[ll].p2;

				if (ptn[p1].y != ptn[p2].y)
					IL = ((double)k - ptn[p2].y) / 
					(ptn[p1].y - ptn[p2].y)
					*nm[idx[p1].org].l + (ptn[p1].y - k) /
					(ptn[p1].y - ptn[p2].y) * nm[idx[p2].org].l;
				else
					IL = nm[idx[p1].org].l;

				IL = fabs(IL);

				p1 = etn[rl].p1;  
				p2 = etn[rl].p2;

				if (ptn[p1].y != ptn[p2].y)
					IR = ((double)k - ptn[p2].y) / 
					(ptn[p1].y - ptn[p2].y) * nm[idx[p1].org].l +
					(ptn[p1].y - k) / 
					(ptn[p1].y - ptn[p2].y) *nm[idx[p2].org].l;
				else
					IR = nm[idx[p2].org].l;

				IR = fabs(IR);

				for (j = floor(minx) + 1; j <= ceil(maxx) - 1; j++)
				{
					dep = (-poly[i].a * 
						(double)j - poly[i].b * (double)k - poly[i].d) 
						/ poly[i].c;

					if (dep>zbuff[j][k])
					{
						if (maxx != minx)
							IP = (maxx - (double)j) / (maxx - minx)*IL +
							((double)j - minx) / (maxx - minx)*IR;
						else 
							IP = IR;

						if (IP>255.0) IP = 255.0;
						else if (IP<0.0) IP = 0.0;

						XSetForeground(display, topwinGc, exact_defs[(int)IP].pixel);
						XDrawPoint(display, toppix, topwinGc, j, k);

						zbuff[j][k] = dep;
					}
				}
			}
		}
	}
}


/*--------------------------------------------------------------------*/
void render_triangles()
{
	double x1[5], y1[5], z1[5];
	double x2[5], y2[5], z2[5];
	double minx, maxx, miny, maxy, xp1, xp2, m, b, xx1, xx2, yy1, yy2;
	double dep;
	int i, j, k, st, end, c;

	for (j = 0; j <= TESTING_HEIGHT; j++)
		for (k = 0; k <= TESTING_WIDTH; k++)
			zbuff[j][k] = -9999999999.0;

	for (i = 0; i <= plycnt; i++)
	{
		st = poly[i].st; 
		end = poly[i].end;

		c = 0;

		miny = (double)TESTING_WIDTH;
		maxy = 0.0;

		for (j = st; j <= end; j++)
			if (etn[j].p1 != -1)
			{
				get_edge2(j, &x1[c], &y1[c], &z1[c], &x2[c], &y2[c], &z2[c]);

				if (y1[c]<miny) miny = y1[c];
				if (y2[c]<miny) miny = y2[c];
				if (y1[c]>maxy) maxy = y1[c];
				if (y2[c]>maxy) maxy = y2[c];

				c++;
			}

		for (k = (int)miny; k <= (int)maxy; k++)
		{
			minx = (double)TESTING_WIDTH;
			maxx = 0.0;

			for (j = 0; j<c; j++)
			{
				if (x1[j]<x2[j]) 
					{ xx1 = x1[j]; xx2 = x2[j]; }
				else 
					{ xx1 = x2[j]; xx2 = x1[j]; }

				if (y1[j]<y2[j]) 
					{ yy1 = y1[j]; yy2 = y2[j]; }
				else 
					{ yy1 = y2[j]; yy2 = y1[j]; }

				if (x2[j] != x1[j])
					m = (y2[j] - y1[j]) / (x2[j] - x1[j]);
				else
					m = 0.000001;

				b = y1[j] - m*x1[j];
				xp1 = ((double)k - b) / m;

				if (xp1<minx && xp1>xx1 &&
					(double)k>yy1 && (double)k<yy2) minx = xp1;

				if (xp1>maxx && xp1<xx2 &&
					(double)k>yy1 && (double)k<yy2) maxx = xp1;
			}

			if (c == 3 && (int)minx>0 && (int)maxx<TESTING_WIDTH)
				for (j = (int)minx; j <= (int)maxx; j++)
				{
					dep = (-poly[i].a*(double)j - poly[i].b*(double)k
						- poly[i].d) / poly[i].c;

					if (dep>zbuff[j][k])
					{
						XSetForeground(display, topwinGc, light(i));
						XDrawPoint(display, toppix, topwinGc, j, k);
						zbuff[j][k] = dep;
					}
				}
		}
	}

}

/*--------------------------------------------------------------------*/
void plane_calc()
{
	int p1, p2, p3;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3, aa, bb, cc, dd;
	int i, st, end, pt;

	for (i = 0; i <= plycnt; i++)
	{
		st = poly[i].st; 
		end = poly[i].end;

		p1 = etn[st].p2;
		p2 = etn[st].p1;
		p3 = etn[st + 1].p1;

		if (p3 == p1 || p3 == p2) 
		{
			p3 = et[st + 1].p2;
		}

		x1 = ptn[p1].x; 
		y1 = ptn[p1].y; 
		z1 = ptn[p1].z;

		x2 = ptn[p2].x; 
		y2 = ptn[p2].y; 
		z2 = ptn[p2].z;

		x3 = ptn[p3].x; 
		y3 = ptn[p3].y; 
		z3 = ptn[p3].z;

		plane(x1, x2, x3, y1, y2, y3, z1, z2, z3, &aa, &bb, &cc, &dd);

		poly[i].a = aa; 
		poly[i].b = bb;
		poly[i].c = cc; 
		poly[i].d = dd;
	}
}

/*--------------------------------------------------------------------*/
void rot_obj(w)
int w;
{
	double x, y, z;
	int i, j, st, end, p1;

	for (i = 0; i <= pntcnt2; i++) 
	{
		roted[i] = 0;
	}

	for (i = obj[w].st; i <= obj[w].end; i++)
	{
		st = poly[i].st; 
		end = poly[i].end;

		for (j = st; j <= end; j++)
		{
			p1 = etn[j].p1;
			if (!roted[p1] && p1 != -1)
			{
				get_pnt(p1, &x, &y, &z);

				rotate(&x, &y, &z, obj[w].rx, obj[w].ry, obj[w].rz);
				translate(&x, &y, &z, obj[w].x, obj[w].y, obj[w].z);

				put_pnt(p1, x, y, z);

				roted[p1] = 1;
			}

			p1 = etn[j].p2;
			if (!roted[p1] && p1 != -1)
			{
				get_pnt(p1, &x, &y, &z);

				rotate(&x, &y, &z, obj[w].rx, obj[w].ry, obj[w].rz);
				translate(&x, &y, &z, obj[w].x, obj[w].y, obj[w].z);

				put_pnt(p1, x, y, z);

				roted[p1] = 1;
			}
		}
	}
}


/*--------------------------------------------------------------------*/
void clip()
{
	int i, code1, code2, ck1, ck2, p1, p2, nxt, tst;
	double u, d1, d2;
	struct point p[1];

	for (i = 0; i <= edgecnt2; i++)
	{
		code1 = 0; 	
		code2 = 0;

		get_edge2(i, &x, &y, &z, &xp, &yp, &zp);

		code1 = getcode(x, y, z); 
		code2 = getcode(xp, yp, zp);

		if (code1 != 0 || code2 != 0)		/* clip! */
			if (code1 & code2)
				etn[i].p1 = -1;		/* outside...ignore! */
			else
			{			/* do clips! */
				nxt = 0;
				u = (Cmin.z - z) / (zp - z);
				if (u >= 0.0 && u <= 1.0)
				{
					p[nxt].x = x + (xp - x)*u;
					p[nxt].y = y + (yp - y)*u;
					p[nxt].z = Cmin.z;

					tst = getcode(p[nxt].x, p[nxt].y, p[nxt].z);

					if (!tst) nxt++;
				}

				u = (Cmax.z - z) / (zp - z);

				if (u >= 0.0 && u <= 1.0)
				{
					p[nxt].x = x + (xp - x)*u;
					p[nxt].y = y + (yp - y)*u;
					p[nxt].z = Cmax.z;

					tst = getcode(p[nxt].x, p[nxt].y, p[nxt].z);

					if (!tst) nxt++;
				}

				u = (Cmax.y - y) / (yp - y);

				if (u >= 0.0 && u <= 1.0)
				{
					p[nxt].x = x + (xp - x)*u;
					p[nxt].y = Cmax.y;
					p[nxt].z = z + (zp - z)*u;

					tst = getcode(p[nxt].x, p[nxt].y, p[nxt].z);

					if (!tst) nxt++;
				}

				u = (Cmin.y - y) / (yp - y);

				if (u >= 0.0 && u <= 1.0)
				{
					p[nxt].x = x + (xp - x)*u;
					p[nxt].y = Cmin.y;
					p[nxt].z = z + (zp - z)*u;

					tst = getcode(p[nxt].x, p[nxt].y, p[nxt].z);

					if (!tst) nxt++;
				}
				u = (Cmin.x - x) / (xp - x);

				if (u >= 0.0 && u <= 1.0)
				{
					p[nxt].x = Cmin.x;
					p[nxt].y = y + (yp - y)*u;
					p[nxt].z = z + (zp - z)*u;

					tst = getcode(p[nxt].x, p[nxt].y, p[nxt].z);

					if (!tst) nxt++;
				}

				u = (Cmax.x - x) / (xp - x);

				if (u >= 0.0 && u <= 1.0)
				{
					p[nxt].x = Cmax.x;
					p[nxt].y = y + (yp - y)*u;
					p[nxt].z = z + (zp - z)*u;

					tst = getcode(p[nxt].x, p[nxt].y, p[nxt].z);

					if (!tst) nxt++;
				}

				ck1 = getcode(p[0].x, p[0].y, p[0].z);
				ck2 = getcode(p[1].x, p[1].y, p[1].z);

				if (nxt == 0 || ck1 && ck2) 
					etn[i].p1 = -1;
				else if (nxt == 1 && !ck1)
					if (code1 && !code2)
					{
						x = p[0].x; 
						y = p[0].y; 
						z = p[0].z;
					}
					else if (code2 && !code1)
					{
						xp = p[0].x; 
						yp = p[0].y; 
						zp = p[0].z;
					}
					else 
						etn[i].p1 = -1;
				else if (nxt == 2 && !ck1 && !ck2)
					etn[i].p1 = -1;
				else 
					etn[i].p1 = -1;
				
				p1 = etn[i].p1; 
				p2 = etn[i].p2;

				if (etn[i].p1 != -1)
				{
					put_pnt(p1, x, y, z); 
					put_pnt(p2, xp, yp, zp);
				}
			}
	}
}


/*--------------------------------------------------------------------*/
int getcode(x, y, z)
double x, y, z;
{
	int code;
	code = 0;

	if (x < Cmin.x) 
	{
		code = code | 1;
	}

	if (x > Cmax.x)
	{
		code = code | 2;
	}

	if (y < Cmin.y)
	{
		code = code | 4;
	}

	if (y > Cmax.y)
	{
		code = code | 8;
	}

	if (z < Cmin.z)
	{
		code = code | 16;
	}

	if (z > Cmax.z)
	{
		code = code | 32;
	}

	return(code);
}


/*--------------------------------------------------------------------*/
void draw_lines_back()
{
	int i, j, st, end;
	for (i = 0; i <= plycnt; i++)
	{
		if (poly[i].c>0.0)
		{
			st = poly[i].st; 
			end = poly[i].end;

			for (j = st; j <= end; j++)
				if (etn[j].p1 != -1)
				{
					get_edge2(j, &x, &y, &z, &xp, &yp, &zp);
					XDrawLine(display, toppix, topwinGc, (int)x, (int)y, (int)xp, (int)yp);
				}
		}
	}
}


/*--------------------------------------------------------------------*/
void draw_lines()
{
	int i;

	for (i = 0; i <= edgecnt2; i++)
	{
		if (etn[i].p1 != -1)
		{
			get_edge2(i, &x, &y, &z, &xp, &yp, &zp);
			XDrawLine(display, toppix, topwinGc, (int)x, (int)y, (int)xp, (int)yp);
		}
	}
}


/*--------------------------------------------------------------------*/
void add_to_list()
{
	int i, j, in, oldin, e1, e2, pc, p1, p2, p3;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3, aa, bb, cc, dd;

	objcnt++;
	obj[objcnt].st = plycnt + 1;

	pc = 0;
	for (i = 0; i <= seg; i++)
	{
		if (i % 1000 == 0) 
		{
			printf("%d ", i / 1000);
		}

		get_line(&a[i], &x, &y, &z, &xp, &yp, &zp);

		in = 1;
		for (j = 0; j <= pntcnt; j++)
		{
			if (x == pt[j].x && y == pt[j].y && z == pt[j].z)
			{
				in = 0;
				oldin = j;
				j = pntcnt + 1;
			}
		}

		if (in)
		{
			pntcnt++;

			pt[pntcnt].x = x;
			pt[pntcnt].y = y;
			pt[pntcnt].z = z;

			e1 = pntcnt;
		}
		else
		{
			e1 = oldin;
		}

		in = 1;
		for (j = 0; j <= pntcnt; j++)
		{
			if (xp == pt[j].x && yp == pt[j].y && zp == pt[j].z)
			{
				in = 0;
				oldin = j;
				j = pntcnt + 1;
			}
		}

		if (in) 
		{
			pntcnt++;

			pt[pntcnt].x = xp; 
			pt[pntcnt].y = yp; 
			pt[pntcnt].z = zp;

			e2 = pntcnt;
		}
		else 
		{
			e2 = oldin;
		}

		edgecnt++;

		et[edgecnt].p1 = e1;
		et[edgecnt].p2 = e2;
		pc++;

		if (pc == 1)
		{
			plycnt++;
			poly[plycnt].st = edgecnt;
		}
		else if (pc == 3)
		{
			poly[plycnt].end = edgecnt; 
			pc = 0;

			p1 = et[poly[plycnt].st].p2;
			p2 = et[poly[plycnt].st].p1;
			p3 = et[poly[plycnt].st + 2].p1;

			if (p3 == p1 || p3 == p2) 
			{
				p3 = et[poly[plycnt].st + 2].p2;
			}
			
			x1 = pt[p1].x; 
			y1 = pt[p1].y; 
			z1 = pt[p1].z;

			x2 = pt[p2].x; 
			y2 = pt[p2].y; 
			z2 = pt[p2].z;

			x3 = pt[p3].x; 
			y3 = pt[p3].y; 
			z3 = pt[p3].z;

			plane(x1, x2, x3, y1, y2, y3, z1, z2, z3, &aa, &bb, &cc, &dd);

			poly[plycnt].a = aa; 
			poly[plycnt].b = bb;
			poly[plycnt].c = cc; 
			poly[plycnt].d = dd;
		}
	}

	obj[objcnt].end = plycnt;
	printf("\n");
}


/*--------------------------------------------------------------------*/
void plane(x1, x2, x3, y1, y2, y3, z1, z2, z3, a, b, c, d)
double x1, x2, x3, y1, y2, y3, z1, z2, z3, *a, *b, *c, *d;
{
	*a = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2);
	*b = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2);
	*c = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2);
	*d = -x1*(y2*z3 - y3*z2) - x2*(y3*z1 - y1*z3) - x3*(y1*z2 - y2*z1);
}


/*--------------------------------------------------------------------*/
double mag(p)
struct point *p;
{ 
	return (sqrt(p->x * p->x + p->y * p->y + p->z * p->z)); 
}


/*--------------------------------------------------------------------*/
void proj_point(PRP, Wmin, Wmax, P)
struct point *PRP, *Wmin, *Wmax, *P;
{
	double a, b;
	struct point n;

	a = -(PRP->x - (Wmin->x + Wmax->x) / 2.0) / PRP->z;	/* shear */
	b = -(PRP->y - (Wmin->y + Wmax->y) / 2.0) / PRP->z;

	n.x = P->x + a * (P->z - PRP->z);
	n.y = P->y + b * (P->z - PRP->z);
	n.z = P->z;
	
	P->x = PRP->x - (n.x - PRP->x) * (PRP->z / (n.z - PRP->z));
	P->y = PRP->y - (n.y - PRP->y) * (PRP->z / (n.z - PRP->z));
}


/*--------------------------------------------------------------------*/
void device_trans_all()
{
	int i;

	for (i = 0; i <= pntcnt2; i++)
	{
		ptn[i].x = TESTING_X + TESTING_WIDTH / 2.0 + ptn[i].x * XSCALE;
		ptn[i].y = TESTING_Y + TESTING_WIDTH / 2.0 + ptn[i].y * YSCALE;
	}
}


/*--------------------------------------------------------------------*/
void proj_trans_all(PRP, Wmin, Wmax)
struct point *PRP, *Wmin, *Wmax;
{
	int i;

	for (i = 0; i <= pntcnt2; i++)
	{
		proj_point(PRP, Wmin, Wmax, &ptn[i]);
	}
}


/*--------------------------------------------------------------------*/
void view_trans_all(VRP, n, u, v)
struct point *VRP, *n, *u, *v;
{
	int i;

	for (i = 0; i <= pntcnt2; i++)
	{
		view_trans(VRP, n, u, v, &ptn[i]);
	}
}


/*--------------------------------------------------------------------*/
void copy_all()
{
	int i, j, st, end;

	pntcnt2 = -1;
	edgecnt2 = -1;

	/* loop over lines and pull out NEW points...set new # of pts */
	for (i = 0; i <= edgecnt; i++)
	{
		get_edge(i, &x, &y, &z, &xp, &yp, &zp);

		edgecnt2++; 
		etn[edgecnt2].p1 = pntcnt2 + 1; 
		etn[edgecnt2].p2 = pntcnt2 + 2;

		pntcnt2++; 
		ptn[pntcnt2].x = x; 
		ptn[pntcnt2].y = y; 
		ptn[pntcnt2].z = z;

		pntcnt2++; 
		ptn[pntcnt2].x = xp; 
		ptn[pntcnt2].y = yp; 
		ptn[pntcnt2].z = zp;
	}
}


/*--------------------------------------------------------------------*/
void view_trans(VRP, n, u, v, P)
struct point *VRP, *n, *u, *v, *P;
{
	struct point p;

	P->x = P->x - VRP->x;		/* TRANSLATE TO VRP */
	P->y = P->y - VRP->y;
	P->z = P->z - VRP->z;

	p.x = P->x * u->x + P->y * u->y + P->z * u->z;
	p.y = P->x * v->x + P->y * v->y + P->z * v->z;
	p.z = P->x * n->x + P->y * n->y + P->z * n->z;

	P->x = p.x;
	P->y = p.y;
	P->z = p.z;
}


/*--------------------------------------------------------------------*/
void find_nuv(VRP, N, V, n, u, v)
struct point *VRP, *N, *V, *n, *u, *v;
{
	double m;

	m = mag(N);
	n->x = N->x / m; 
	n->y = N->y / m; 
	n->z = N->z / m;
	cross(V, N, u);

	m = mag(u);
	u->x = u->x / m; 
	u->y = u->y / m; 
	u->z = u->z / m;
	cross(n, u, v);
}


/*--------------------------------------------------------------------*/
void cross(p1, p2, pc)
struct point *p1, *p2, *pc;
{
	pc->x = p1->y * p2->z - p1->z * p2->y;
	pc->y = p1->z * p2->x - p1->x * p2->z;
	pc->z = p1->x * p2->y - p1->y * p2->x;
}


/*--------------------------------------------------------------------*/
void translate(x, y, z, ax, ay, az)
double *x, *y, *z;
double ax, ay, az;
{
	*x = *x + ax;
	*y = *y + ay;
	*z = *z + az;
}


/*--------------------------------------------------------------------*/
void rotate(x, y, z, ax, ay, az)
double *x, *y, *z;
double ax, ay, az;
{
	double xp, yp, zp;

	xp = *x * cos(az) - *y * sin(az);	/* Y AXIS */
	yp = *x * sin(az) + *y * cos(az);
	*x = xp;	
	*y = yp;

	yp = *y * cos(ax) - *z * sin(ax);	/* Z AXIS */
	zp = *y * sin(ax) + *z * cos(ax);
	*y = yp;	
	*z = zp;
	
	zp = *z * cos(ay) - *x * sin(ay);	/*  AXIS */
	xp = *z * sin(ay) + *x * cos(ay);
	*z = zp;	
	*x = xp;
}


/*--------------------------------------------------------------------*/
double pw(x, y)
double x, y;
{
	double res;
	int neg;

	if (x < 0)
	{
		neg = 1;
	}
	else
	{
		neg = 0;
	}

	if (x == 0.0)
	{
		x = 0.000000000001;
	}

	res = exp(y*log(fabs(x)));

	if (neg)
	{
		return(-res);
	}
	else
	{
		return(res);
	}
}


/*--------------------------------------------------------------------*/
double to_rad(d)
double d;
{ 
	return(2 * 3.14159 * (d / 180.0)); 
}


/*---------------------------------------------------------------------------*/
void get_point(i, x1, y1, z1)
int i;
double *x1, *y1, *z1;
{
	*x1 = pt[i].x;
	*y1 = pt[i].y;
	*z1 = pt[i].z;
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void get_pnt(i, x, y, z)
int i;
double *x, *y, *z;
{
	*x = ptn[i].x;
	*y = ptn[i].y;
	*z = ptn[i].z;
}


/*---------------------------------------------------------------------------*/
void put_pnt(i, x, y, z)
int i;
double x, y, z;
{
	ptn[i].x = x;
	ptn[i].y = y;
	ptn[i].z = z;
}


/*---------------------------------------------------------------------------*/
void get_edge2(i, x1, y1, z1, x2, y2, z2)
int i;
double *x1, *y1, *z1, *x2, *y2, *z2;
{
	int v1, v2;

	v1 = etn[i].p1;  v2 = etn[i].p2;

	*x1 = ptn[v1].x; *x2 = ptn[v2].x;
	*y1 = ptn[v1].y; *y2 = ptn[v2].y;
	*z1 = ptn[v1].z; *z2 = ptn[v2].z;
}


/*---------------------------------------------------------------------------*/
void get_edge(i, x1, y1, z1, x2, y2, z2)
int i;
double *x1, *y1, *z1, *x2, *y2, *z2;
{
	int v1, v2;

	v1 = et[i].p1;  
	v2 = et[i].p2;

	*x1 = pt[v1].x; 
	*x2 = pt[v2].x;
	*y1 = pt[v1].y; 
	*y2 = pt[v2].y;
	*z1 = pt[v1].z; 
	*z2 = pt[v2].z;
}


/*---------------------------------------------------------------------------*/
void get_line(p, x1, y1, z1, x2, y2, z2)
struct line *p;
double *x1, *y1, *z1, *x2, *y2, *z2;
{
	*x1 = p->x1; 
	*x2 = p->x2; 
	*y1 = p->y1; 
	*y2 = p->y2; 
	*z1 = p->z1; 
	*z2 = p->z2; 
}


/*---------------------------------------------------------------------------*/
void set_line(p, x1, y1, z1, x2, y2, z2)
struct line *p;
double x1, y1, z1, x2, y2, z2;
{ 
	p->x1 = x1; 
	p->x2 = x2; 
	p->y1 = y1; 
	p->y2 = y2; 
	p->z1 = z1; 
	p->z2 = z2; 
}


/*---------------------------------------------------------------------------*/
void set_point(p, x, y, z)
struct point *p;
double x, y, z;
{ 
	p->x = x; 
	p->y = y; 
	p->z = z; 
}


/*---------------------------------------------------------------------------*/
/*	
case 1: XDrawLine(display, topwin, topwinGc, 20, 20, 100, 100);
case 3: XDrawRectangle(display, topwin, topwinGc, 100, 150, 20, 40);
case 4: XFillRectangle(display, topwin, topwinGc, 300, 300, 50, 10);
case 5: XDrawPoint(display, topwin, topwinGc, 400, 400);
*/
/*---------------------------------------------------------------------------*/

void load_object(name)
char name[];
{
	FILE *fp;
	double xf, yf, zf, scale;
	int num, i, j;

	seg = -1;   
	j = 0;
	
	printf("enter a obj file >");
	gets(name);

	fp = fopen(name, "r");
	if (fp == NULL)
	{
		exit(-1);
	}

	scale = 0.00000001;

	while (!feof(fp))
	{
		fscanf(fp, "%d\n", &num);
		for (i = 1; i <= num; i++)
		{
			j++;
			fscanf(fp, "%lf%lf%lf\n", &x, &y, &z);

			if (fabs(x)>scale) scale = fabs(x);
			if (fabs(y)>scale) scale = fabs(y);
			if (fabs(z)>scale) scale = fabs(z);

			if (i == 1)
			{
				xf = x;
				yf = y;
				zf = z;
			}
			else
			{
				seg++;
				set_line(&a[seg], x, y, z, xp, yp, zp);
			}

			xp = x;
			yp = y;
			zp = z;
		}

		seg++;
		set_line(&a[seg], xp, yp, zp, xf, yf, zf);
		
		if (j % 100 == 0)
		{
			printf("%d ", j / 100);
		}
	}
	printf("\n");

	scale = SCALEFAC / scale;
	for (i = 0; i <= seg; i++)
	{
		a[i].x1 *= scale; a[i].x2 *= scale;
		a[i].y1 *= scale; a[i].y2 *= scale;
		a[i].z1 *= scale; a[i].z2 *= scale;
	}

	printf("***BUILD***\n");
	add_to_list();
	fclose(fp);
}

/*---------------------------------------------------------------------------*/
void make_elipse()
{
	printf("ELLIPSOID GENERATOR\n---------------------------\n");
	printf("Enter s1 >");
	scanf("%lf", &s1);
	printf("Enter s2 >");
	scanf("%lf", &s2);
	printf("Enter rx >");
	scanf("%lf", &rx);
	printf("Enter ry >");
	scanf("%lf", &ry);
	printf("Enter rz >");
	scanf("%lf", &rz);
	printf("Enter Number of ELEMENTS >");
	scanf("%d", &DIV);

	printf("** BUILDING **\n\n");
	seg = -1;

	for (g = -PI / 2; g <= PI / 2 + .001; g = g + PI / DIV)
	{
		t = -PI;

		x = (rx * pw(cos(g), s1) * pw(cos(t), s2));
		y = (ry * pw(cos(g), s1) * pw(sin(t), s2));
		z = (rz * pw(sin(g), s1));

		sv = -1;

		for (t = -PI; t <= PI; t = t + PI / DIV)
		{
			xp = x; yp = y; zp = z;

			x = (rx * pw(cos(g), s1) * pw(cos(t), s2));
			y = (ry * pw(cos(g), s1) * pw(sin(t), s2));
			z = (rz * pw(sin(g), s1));

			seg++;
			set_line(&a[seg], xp, yp, zp, x, y, z);

			sv++;
			if (g != -PI / 2)
			{
				seg++;
				set_line(&a[seg], x, y, z,
					bo[sv].x2, bo[sv].y2, bo[sv].z2);
			}

			set_line(&b[sv], xp, yp, zp, x, y, z);
		}

		seg++;
		set_line(&a[seg], b[0].x2, b[0].y2, b[0].z2,
			b[sv].x2, b[sv].y2, b[sv].z2);

		for (i = 0; i <= sv; i++)
		{
			set_line(&bo[i], b[i].x1, b[i].y1, b[i].z1,
				b[i].x2, b[i].y2, b[i].z2);
		}
	}

	add_to_list();
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void init()
{
	if (!(display = XOpenDisplay("")))
	{
		fprintf(stderr, "Unable to open display...exiting !\n");
		exit(-1);
	}

	screen_num = XDefaultScreen(display);
	default_visual = DefaultVisual(display, screen_num);
	default_cmap = XDefaultColormap(display, screen_num);

	NumColors = XDisplayCells(display, screen_num);
	NumColors = 64;
	printf("COLORS! = %ld\n", NumColors);

	/* set colors - Black and White */
	fgColor = WhitePixel(display, screen_num);
	bgColor = BlackPixel(display, screen_num);
	bdColor = WhitePixel(display, screen_num);

	for (i = 0; i<NumColors; i++)
	{
		color.pixel = i;
		XQueryColor(display, default_cmap, &color);

		color.red = i * 255 * 4;
		color.green = i * 255 * 4;
		color.blue = i * 255 * 4;

		if (!XAllocColor(display, default_cmap, &color))
		{
			printf("cant alloc color\n");
		}
		exact_defs[i].pixel = color.pixel;
	}

	/* load the font for the window - exits if not found */
	if ((fontSt = XLoadQueryFont(display, TESTING_DEFAULTFONT)) == NULL)
	{
		fprintf(stderr, "Unable to load font...exiting !\n");
		exit(-1);
	}

	/* create the window with some attr. */
	topwin = XCreateSimpleWindow(
		display,
		DefaultRootWindow(display),
		TESTING_X, TESTING_Y,
		TESTING_WIDTH, TESTING_HEIGHT,
		TESTING_BDWIDTH, bdColor, bgColor
	);

	/* create the pixmap */
	toppix = XCreatePixmap(
		display,
		DefaultRootWindow(display),
		TESTING_WIDTH, TESTING_HEIGHT,
		DefaultDepth(display, screen_num)
	);

	/* set an event of interest for windows */
	XSelectInput(display, topwin, KeyPressMask | ButtonPressMask);

	/* create a graphics context */
	topwinGc = XCreateGC(display, topwin, 0, 0);

	XSetFont(display, topwinGc, fontSt->fid);
	XSetForeground(display, topwinGc, fgColor);

	XMapWindow(display, topwin);

	XSetForeground(display, topwinGc, bgColor);

	XFillRectangle(display, toppix, topwinGc, 0, 0,
		TESTING_WIDTH, TESTING_HEIGHT);

	XSetForeground(display, topwinGc, fgColor);
}


/*--------------------------------------------------------------------*/
void close_all()
{
	XFreeGC(display, topwinGc);
	XFreeFont(display, fontSt);

	XDestroyWindow(display, topwin);
	XCloseDisplay(display);
	exit(0);
}

