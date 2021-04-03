#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NX 200
#define NY 200
#define DT 0.25
#define RE 15000
#define DATA_SIZE (sizeof(double)*NX*NY)
#define PIX_SIZE (3 * NY * NX)

#define C(array, x, y) array[NY*(x) + (y)]
#define R(array, x, y) array[NY*(x+1) + (y)]
#define L(array, x, y) array[NY*(x-1) + (y)]
#define T(array, x, y) array[NY*(x) + (y+1)]
#define B(array, x, y) array[NY*(x) + (y-1)]
#define R2(array, x, y) array[NY*(x+2) + (y)]
#define L2(array, x, y) array[NY*(x-2) + (y)]
#define T2(array, x, y) array[NY*(x) + (y+2)]
#define B2(array, x, y) array[NY*(x) + (y-2)]

double *p_cur;
double *u_cur;
double *v_cur;
double *p_pre;
double *u_pre;
double *v_pre;
unsigned char *pix;

double updatePressure(int x, int y);
void updateVelocity(int x, int y);
double schemeX(double u, double *f, int x, int y);
double schemeY(double u, double *f, int x, int y);
void setWall();
void setObs();
void genWind(int t);
void toHue(unsigned char *pix, double value, double power);

#pragma pack(1)
typedef struct {
	char type[2];
	unsigned int fileSize;
	unsigned short reserved1;
	unsigned short reserved2;
	unsigned int ofs;

	unsigned int chkSize;
	int width;
	int height;
	unsigned short planes;
	unsigned short bits;
	unsigned int comp;
	unsigned int dataSize;
	unsigned int dpmX;
	unsigned int dpmY;
	unsigned int color;
	unsigned int important;
} BMP_HEADER;
#pragma pack()

int main() {
    FILE *fp = NULL;
    char fname[256] = { 0 };
    BMP_HEADER hd = { 0 };

    double resi = 0.0;
    int t = 0, q = 0;
    int x, y;

    hd.chkSize = 40;
    hd.bits = 24;
    hd.type[0] = 'B';
    hd.type[1] = 'M';
    hd.planes = 1;
    hd.width = NX;
    hd.height = NY;
    hd.ofs = 54;
    hd.dataSize = PIX_SIZE;
    hd.fileSize = hd.dataSize + hd.ofs;

    p_cur = (double*)malloc(DATA_SIZE);
    u_cur = (double*)malloc(DATA_SIZE);
    v_cur = (double*)malloc(DATA_SIZE);
    p_pre = (double*)malloc(DATA_SIZE);
    u_pre = (double*)malloc(DATA_SIZE);
    v_pre = (double*)malloc(DATA_SIZE);
    memset(p_cur, 0, DATA_SIZE);
    memset(u_cur, 0, DATA_SIZE);
    memset(v_cur, 0, DATA_SIZE);
    memset(p_pre, 0, DATA_SIZE);
    memset(u_pre, 0, DATA_SIZE);
    memset(v_pre, 0, DATA_SIZE);

    pix = (unsigned char*)malloc(PIX_SIZE);
    unsigned char *ppix;
    double theta;
    double power;
    for (t = 0; q < 128; ++t) {
        if (t % 20 == 0) {
            sprintf_s(fname, sizeof(fname), "C:\\Users\\owner\\source\\repos\\CFD\\data\\%d.bmp", t / 20);
            fopen_s(&fp, fname, "wb");
            fwrite(&hd, sizeof(hd), 1, fp);
            memset(pix, 0, PIX_SIZE);
        }

        setWall();
        setObs();
        genWind(t);

        for (q = 0; q < 128; ++q) {
            resi = 0.0;
            memcpy(p_pre, p_cur, DATA_SIZE);
            for (x = 2; x < NX - 2; ++x) {
                for (y = 2; y < NY - 2; ++y) {
                    resi += updatePressure(x, y);
                }
            }
            if (resi <= 0.25) {
                break;
            }
        }

        memcpy(u_pre, u_cur, DATA_SIZE);
        memcpy(v_pre, v_cur, DATA_SIZE);
        for (x = 2; x < NX - 2; ++x) {
            ppix = pix + 3 * NY * x + 6;
            for (y = 2; y < NY - 2; ++y) {
                updateVelocity(x, y);
                power = sqrt(C(u_cur, x, y)*C(u_cur, x, y) + C(v_cur, x, y)*C(v_cur, x, y));
                theta = atan2(C(v_cur, x, y), C(u_cur, x, y));
                if (theta < 0) { theta += 8 * atan(1); }
                theta /= 8 * atan(1);
                toHue(ppix, theta, power);
                ppix += 3;
            }
        }

        printf("%05d %03d\b\b\b\b\b\b\b\b\b", t, q);

        if (t % 20 == 0) {
            fwrite(pix, PIX_SIZE, 1, fp);
            fclose(fp);
        }
    }

    free(p_cur);
    free(u_cur);
    free(v_cur);
    free(p_pre);
    free(u_pre);
    free(v_pre);
    free(pix);

    return 0;
}

double updatePressure(int x, int y) {
    double ux = (R(u_cur, x, y) - L(u_cur, x, y)) / 2.0;
    double uy = (T(u_cur, x, y) - B(u_cur, x, y)) / 2.0;
    double vx = (R(v_cur, x, y) - L(v_cur, x, y)) / 2.0;
    double vy = (T(v_cur, x, y) - B(v_cur, x, y)) / 2.0;
    double dp = ((R(p_pre, x, y) + L(p_pre, x, y) + T(p_pre, x, y) + B(p_pre, x, y))
        + (ux*ux + vy*vy + 2.0*uy*vx - (ux + vy) / DT)
    ) / 4.0 - C(p_pre, x, y);

    C(p_cur, x, y) = C(p_pre, x, y) + dp * 0.25;

    return dp * dp;
}

void updateVelocity(int x, int y) {
	C(u_cur, x, y) = C(u_pre, x, y) + DT * (
        (L(u_pre, x, y) - 4.0 * C(u_pre, x, y) + R(u_pre, x, y) + B(u_pre, x, y) + T(u_pre, x, y)) / RE
        - (schemeX(C(u_pre, x, y), u_pre, x, y) + schemeY(C(v_pre, x, y), u_pre, x, y))
        + (L(p_cur, x, y) - R(p_cur, x, y)) / 2.0
    );
    C(v_cur, x, y) = C(v_pre, x, y) + DT * (
        (L(v_pre, x, y) - 4.0 * C(v_pre, x, y) + R(v_pre, x, y) + B(v_pre, x, y) + T(v_pre, x, y)) / RE
        - (schemeX(C(u_pre, x, y), v_pre, x, y) + schemeY(C(v_pre, x, y), v_pre, x, y))
        + (B(p_cur, x, y) - T(p_cur, x, y)) / 2.0
    );
}

double schemeX(double u, double *f, int x, int y) {
    return u * (L2(f, x, y) - 8.0 * (L(f, x, y) - R(f, x, y)) - R2(f, x, y)) / 12.0
        + fabs(u) * (L2(f, x, y) - 4.0 * L(f, x, y) + 6.0 * C(f, x, y) - 4.0 * R(f, x, y) + R2(f, x, y)) / 12.0
    ;
}

double schemeY(double u, double *f, int x, int y) {
    return u * (B2(f, x, y) - 8.0 * (B(f, x, y) - T(f, x, y)) - T2(f, x, y)) / 12.0
        + fabs(u) * (B2(f, x, y) - 4.0 * B(f, x, y) + 6.0 * C(f, x, y) - 4.0 * T(f, x, y) + T2(f, x, y)) / 12.0
    ;
}

void setWall() {
	int x, y;

	for (x = 0; x < NX; ++x) {
		C(p_cur, x, 0) = 0.0;
		C(p_cur, x, 1) = 0.0;
		C(p_cur, x, NY - 2) = 0.0;
		C(p_cur, x, NY - 1) = 0.0;

		C(u_cur, x, 0) = 0.0;// P(curU, x, 2);
		C(u_cur, x, 1) = 0.0;// P(curU, x, 2);
		C(u_cur, x, NY - 2) = 0.0;// P(curU, x, NY - 3);
		C(u_cur, x, NY - 1) = 0.0;// P(curU, x, NY - 3);

		C(v_cur, x, 0) = 0.0;// P(curV, x, 2);
		C(v_cur, x, 1) = 0.0;// P(curV, x, 2);
		C(v_cur, x, NY - 2) = 0.0;// P(curV, x, NY - 3);
		C(v_cur, x, NY - 1) = 0.0;// P(curV, x, NY - 3);
	}

	for (y = 0; y < NY; ++y) {
		C(p_cur, 0, y) = 0.0;
		C(p_cur, 1, y) = 0.0;
		C(p_cur, NX - 2, y) = 0.0;
		C(p_cur, NX - 1, y) = 0.0;

		C(u_cur, 0, y) = 0.0;// P(curU, 2, y);
		C(u_cur, 1, y) = 0.0;// P(curU, 2, y);
		C(u_cur, NX - 2, y) = 0.0;// P(curU, NX - 3, y);
		C(u_cur, NX - 1, y) = 0.0;// P(curU, NX - 3, y);

		C(v_cur, 0, y) = 0.0;// P(curV, 2, y);
		C(v_cur, 1, y) = 0.0;// P(curV, 2, y);
		C(v_cur, NX - 2, y) = 0.0;// P(curV, NX - 3, y);
		C(v_cur, NX - 1, y) = 0.0;// P(curV, NX - 3, y);
	}
}

void setObs() {
	int x, y;
	double px, py;
	double r;

	for (x = 0; x < NX; ++x) {
		px = 2.0*x / NX - 1.3;
		for (y = 0; y < NY; ++y) {
			py = 2.0*y / NY - 1.0;
			r = sqrt(px*px + py * py);

			if (0.1 <= r && r <= 0.2) {
				C(p_cur, x, y) = 0.0;
				C(u_cur, x, y) = 0.0;
				C(v_cur, x, y) = 0.0;
			}
		}
	}
}

void genWind(int t) {
    int y;
    double u, v;

    if ((int)(t*DT) % 4000 < 2000) {
        u = 0.5;
        v = 0.0;// 0.25;
    } else {
        u = 0.5;
        v = 0.0;// -0.25;
    }

    for (y = 7 * NY / 16; y < 9 * NY / 16; ++y) {
        C(p_cur, NX / 8, y) = 1.0;
        C(u_cur, NX / 8, y) = u;
        C(v_cur, NX / 8, y) = v;
    }
}

void toHue(unsigned char *pix, double value, double power) {
    value *= 1535;

    if (power < 0) { power = 0; }
    if (1 < power) { power = 1; }

    // Blue
    if (value < 0) {
        pix[0] = (unsigned char)(power * 255);
        pix[1] = 0;
        pix[2] = 0;
        return;
    }
    // Blue - SkyBlue
    if (value < 256) {
        pix[0] = (unsigned char)(power * 255);
        pix[1] = (unsigned char)(power * value);
        pix[2] = 0;
        return;
    }
    // SkyBlue - Green
    if (value < 512) {
        pix[0] = (unsigned char)(power * (511 - value));
        pix[1] = (unsigned char)(power * 255);
        pix[2] = 0;
        return;
    }
    // Green - Yellow
    if (value < 768) {
        pix[0] = 0;
        pix[1] = (unsigned char)(power * 255);
        pix[2] = (unsigned char)(power * (value - 512));
        return;
    }
    // Yellow - Red
    if (value < 1024) {
        pix[0] = 0;
        pix[1] = (unsigned char)(power * (1023 - value));
        pix[2] = (unsigned char)(power * 255);
        return;
    }
    // Red - Purple
    if (value < 1280) {
        pix[0] = (unsigned char)(power * (value - 1024));
        pix[1] = 0;
        pix[2] = (unsigned char)(power * 255);
        return;
    }
    // Purple - Blue
    if (value < 1536) {
        pix[0] = (unsigned char)(power * 255);
        pix[1] = 0;
        pix[2] = (unsigned char)(power * (1535 - value));
        return;
    }

    pix[0] = (unsigned char)(power * 255);
    pix[1] = (unsigned char)(power * 255);
    pix[2] = (unsigned char)(power * 255);
}
