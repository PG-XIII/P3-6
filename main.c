#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
//#include <GL/glut.h>
#include <GLUT/glut.h> // Essa é a versão do glut.h no mac
#include <errno.h>

//rotinas auxiliares
//funcoes que possuem ret como parametro usam ele como destino do resultado da funcao.
float prod_interno(float vec1[], float vec2[]);
void normalizar(float vec[], float ret[]);
void prod_vetorial(float vec1[], float vec2[], float ret[]);
void proj_vetores(float vec1[], float vec2[], float ret[]);
void mul_escalar(float vec[], float k, float ret[]);
void sub_vet(float vec1[], float vec2[], float ret[]);
void sum_vet(float vec1[], float vec2[], float ret[]);
void projetar_triangulo(int *triangulo, float **projecao);
void scanline(float **projecao, int** ret);

void carregar_camera();
void carregar_iluminacao();
void carregar_objetos();
void normalizar_triangulos();
void normalizar_vertices();

void mudanca_base_scc(float vec[], float ret[]);

//variaveis da camera
float C[3];
float N[3];
float V[3];
float U[3];
float d;
float hx;
float hy;

//variaveis de iluminacao
float Pl[3];
float ka;
float Ia[3];
float kd;
float Od[3];
float ks;
float Il[3];
float n;

//variaveis de objeto
int num_pontos;
int num_triangulos;
float** pontos;
int** triangulos;
float** normais_vertices;
float** normais_triangulos;

// Variáveis do z-buffer
int** ds;
int** cores;

// Variáveis da interface gráfica
int width = 640;
int height = 320;

int main(int argc, char **argv)
{
    carregar_objetos();
    normalizar_triangulos();
    normalizar_vertices();
    int i;
    for(i = 0; i < num_triangulos; i++)
    {
        printf("%f %f %f\n", normais_vertices[i][0], normais_vertices[i][1], normais_vertices[i][2]);
    }
    
    return 0;
}

void carregar_objetos()
{
    FILE *fp;
    fp = fopen("objeto.byu","r");
    if (fp == NULL) {
        printf ("Erro de leitura em arquivo, errno = %d\n", errno);
    }

    fscanf(fp," %d %d", &num_pontos, &num_triangulos);

    int i;
    pontos =  (float**) calloc(num_pontos, sizeof(float));
    for(i = 0; i < num_pontos; i++)
    {
        pontos[i] = (float*) calloc(3, sizeof(float));
        fscanf(fp," %f %f %f", &pontos[i][0], &pontos[i][1], &pontos[i][2]);
    }

    triangulos = (int**) calloc(num_triangulos, sizeof(int));
    for(i = 0; i < num_triangulos; i++)
    {
        triangulos[i] = (int*) calloc(3, sizeof(int));
        fscanf(fp," %d %d %d", &triangulos[i][0], &triangulos[i][1], &triangulos[i][2]);
    }

    fclose (fp);
}

void normalizar_triangulos()
{
    int i;
    float aux1[3];
    float aux2[3];
    float aux3[3];
    
    normais_triangulos = (float**) calloc(num_triangulos, sizeof(float));
    for(i = 0; i < num_triangulos; i++)
    {
        sub_vet(pontos[triangulos[i][1]-1], pontos[triangulos[i][0]-1], aux1);
        sub_vet(pontos[triangulos[i][2]-1], pontos[triangulos[i][0]-1], aux2);
        prod_vetorial(aux1, aux2, aux3);
        normalizar(aux3, aux3);
        normais_triangulos[i] = (float*) calloc(3, sizeof(float));
        normais_triangulos[i][0] = aux3[0];
        normais_triangulos[i][1] = aux3[1];
        normais_triangulos[i][2] = aux3[2];
    }
}

void normalizar_vertices()
{
    int i,j,k;
    float inc = 0;
    float aux1[3] = {0,0,0};

    normais_vertices = (float**) calloc(num_pontos, sizeof(float));
    for(i = 0; i < num_pontos; i++)
    {
        for(j = 0; j < num_triangulos; j++)
        {
            for(k = 0; k < 3; k++)
            {
                if(triangulos[j][k] == (i + 1))
                {
                    sum_vet(aux1, normais_triangulos[j], aux1);
                    inc++;
                }
            }
        }

        mul_escalar(aux1, 1/inc, aux1);
        normalizar(aux1, aux1);
        normais_vertices[i] = (float*) calloc(3, sizeof(float));
        normais_vertices[i][0] = aux1[0];
        normais_vertices[i][1] = aux1[1];
        normais_vertices[i][2] = aux1[2];
        inc = 0;
        aux1[0] = 0;
        aux1[1] = 0;
        aux1[2] = 0;
    }
}

void carregar_camera()
{
    FILE *fp;
    fp = fopen("camera.cfg","r");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
    }
    fscanf(fp," %f %f %f", &C[0], &C[1], &C[2]);
    fscanf(fp," %f %f %f", &N[0], &N[1], &N[2]);
    fscanf(fp," %f %f %f", &V[0], &V[1], &V[2]);
    fscanf(fp," %f %f %f", &d, &hx, &hy);

    fclose (fp);

    //Ortogonalizar V
    float aux1[3];
    float aux2[3];
    normalizar(N,N);
    proj_vetores(V, N, aux1);
    sub_vet(V, aux1, aux2);
    V[0] = aux2[0];
    V[1] = aux2[1];
    V[2] = aux2[2];
    normalizar(V,V);
    //Encontrar U
    prod_vetorial(V, N, U);
}

void carregar_iluminacao()
{
    FILE *fp;
    fp = fopen("iluminacao.txt","r");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
    }
    fscanf(fp," %f %f %f", &Pl[0], &Pl[1], &Pl[2]);
    fscanf(fp," %f", &ka);
    fscanf(fp," %f %f %f", &Ia[0], &Ia[1], &Ia[2]);
    fscanf(fp," %f", &kd);
    fscanf(fp," %f %f %f", &Od[0], &Od[1], &Od[2]);
    fscanf(fp," %f", &ks);
    fscanf(fp," %f %f %f", &Il[0], &Il[1], &Il[2]);
    fscanf(fp," %f", &n);

    fclose (fp);
}

void projetar_triangulo(int *triangulo, float **projecao) {
    int i;
    for (i = 0; i < 3; i++) {
        float* p = pontos[triangulo[i]];
        float pbarra[3]; // Ponto em coordenadas de câmera

        mudanca_base_scc(p, pbarra);

        projecao[i][0] = ((pbarra[0] * d / pbarra[2] / hx) + 1)*width/2; // Esse valor precisa apenas ser multiplicado pela [largura/altura] da tela em que será apresentado e arredondado
        projecao[i][1] = (1 - (pbarra[1] * d / pbarra[2] / hy))*height/2;
    }
}

void scanline(float** projecao, int** ret) {
    float* top = projecao[0];
    float* middle = projecao[1];
    float* bottom = projecao[2];

    // Calcular o primeiro pedaço dos xmin e xmax
    float a = (bottom[1]-top[1])/(bottom[0]-top[0]);
    float b = (middle[1]-top[1])/(middle[0]-top[0]);
    int i, j;
    for (i = floor(top[1]), j = 0; i < floor(middle[1]); i++, j++) {
        int x1 = floor(top[0] + j/a);
        int x2 = floor(top[0] + j/b);

        if (x1 > x2) {
            int tmp = x1;
            x1 = x2;
            x2 = tmp;
        }

        ret[j][0] = x1;
        ret[j][1] = x2;
    }

    // Calcular o segundo pedaço
    b = (bottom[1]-middle[1])/(bottom[0]-middle[0]);
    for (; i < bottom[1]; i++, j++) {
        int x1 = floor(top[0] + j/a);
        int x2 = floor(middle[0] + (j-i)/b);

        if (x1 > x2) {
            int tmp = x1;
            x1 = x2;
            x2 = tmp;
        }

        ret[j][0] = x1;
        ret[j][1] = x2;
    }
}

/* Funções do z-buffer */

void preencher_z_buffer() {
    int i;
    for (i = 0; i < num_triangulos; i++) {
        float **projecao = (float**)malloc(3*sizeof(float*));
        int i;
        for (i = 0; i < 3; i++) {
            projecao[i] = (float*)malloc(2*sizeof(float));
        }
        projetar_triangulo(triangulos[i], projecao);
        
        // Ordenar os pontos pela coordenada y
        float max = -999999999;
        int max_i = 0;
        float min = 999999999;
        int min_i = 0;
        for (i = 0; i < 3; i++) {
            if (projecao[i][1] > max) {
                max_i = i;
                max = projecao[i][1];
            } 
       
            if (projecao[i][1] < min) {
                min_i = i;
                min = projecao[i][1];
            }
        }

        float* top = projecao[max_i];
        float* bottom = projecao[min_i];
        float* middle = projecao[0 + 1 + 2 - max_i - min_i];

        projecao[0] = top;
        projecao[1] = middle;
        projecao[2] = bottom;

        int n_linhas = floor(max) - floor(min);
        int **xminmax = (int**)malloc(n_linhas * sizeof(int*));
        for (i = 0; i < n_linhas; i++) {
            xminmax[i] = (int*)malloc(2*sizeof(int));
        }

        scanline(projecao, xminmax);

        int j;
        for (i = 0; i < n_linhas; i++) {
            for (j = xminmax[i][0]; j < xminmax[i][1]; j++) {
                
            }
        }
    }
}

/* Funções Algébricas */

float prod_interno(float vec1[], float vec2[])
{
    return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]);
}

void normalizar(float vec[], float ret[])
{
    float mod = sqrt(prod_interno(vec,vec));

    ret[0] = vec[0]/mod;
    ret[1] = vec[1]/mod;
    ret[2] = vec[2]/mod;
}

void prod_vetorial(float vec1[], float vec2[], float ret[])
{
    ret[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    ret[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    ret[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

void proj_vetores(float vec1[], float vec2[], float ret[])
{
    float prop = prod_interno(vec1,vec2)/prod_interno(vec2,vec2);

    ret[0] = prop * vec2[0];
    ret[1] = prop * vec2[1];
    ret[2] = prop * vec2[2];
}

void mul_escalar(float vec[], float k, float ret[])
{
    ret[0] = vec[0] * k;
    ret[1] = vec[1] * k;
    ret[2] = vec[2] * k;
}

void sub_vet(float vec1[], float vec2[], float ret[])
{
    ret[0] = vec1[0] - vec2[0];
    ret[1] = vec1[1] - vec2[1];
    ret[2] = vec1[2] - vec2[2];
}

void sum_vet(float vec1[], float vec2[], float ret[])
{
    ret[0] = vec1[0] + vec2[0];
    ret[1] = vec1[1] + vec2[1];
    ret[2] = vec1[2] + vec2[2];
}

void mudanca_base_scc(float vec[], float ret[])
{
    ret[0] = U[0]*(vec[0] - C[0]) + U[1]*(vec[1] - C[1]) + U[2]*(vec[2] - C[2]);
    ret[1] = V[0]*(vec[0] - C[0]) + V[1]*(vec[1] - C[1]) + V[2]*(vec[2] - C[2]);
    ret[2] = N[0]*(vec[0] - C[0]) + N[1]*(vec[1] - C[1]) + N[2]*(vec[2] - C[2]);
}
