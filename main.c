#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#ifdef __APPLE__
    // Essa é a versão da glut pra macOS
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif
#include <errno.h>
#include <unistd.h>

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
void coordenadas_baricentricas(int* ponto, float** triangulo, float* coordenadas);
void resolver_sistema(float** matriz, int n, int m, float* resultado);
void escalonar(float** matriz, int n, int m);

void draw();
void carregar_camera();
void carregar_iluminacao();
void iluminar(float V[], float N[], float L[], float ret[]);
void carregar_objetos();
void normalizar_triangulos();
void coord_mundo_para_scc();
void normalizar_vertices();
void projetar_pontos();
void init_z_buffer();
void preencher_z_buffer();
void scanline(int **projecao, int** ret);


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
int** pontos_projetados;

// Variáveis do z-buffer
float** z_buffer_d;
float*** z_buffer_cor;

// Resolucao da interface gráfica
#define width 500
#define height 500

int main(int argc, char **argv)
{
    printf("Carregando camera...");
    carregar_camera();
    printf("OK\nCarregando objetos...");
    carregar_objetos();
    printf("OK\nCarregar iluminação...");
    carregar_iluminacao();
    printf("OK\nMudando as coordenadas de mundo para coordenadas de câmera...");
    coord_mundo_para_scc();
    printf("OK\nEncontrando pontos em coordenadas de tela...");
    projetar_pontos();
    printf("OK\nNormalizando triangulos...");
    normalizar_triangulos();
    printf("OK\nNormalizando vertices...");
    normalizar_vertices();
    printf("OK\nInicializando o z-buffer...");
    init_z_buffer();
    printf("OK\nPreenchendo o z-buffer...\n");
    preencher_z_buffer();
    printf("                         OK\n");

    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_RGB);
    glutInitWindowSize(width,height);
    glutInitWindowPosition(100,100);
    glutCreateWindow("PG-13");
    glClearColor(0.0,0.0,0.0,0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    gluOrtho2D(0.0, width, 0.0, height);

/*    glBegin(GL_POINTS);
        glColor3f(1,1,1);
        int i;
        for(i=0; i < num_pontos; i++)
        {
            glVertex2i(pontos_projetados[i][0], pontos_projetados[i][1]);
        }
    glEnd();*/
    glFlush();
    glutDisplayFunc(draw);
    glutMainLoop();

    free(pontos);
    free(triangulos);
    free(normais_vertices);
    free(normalizar_triangulos);
    free(pontos_projetados);
    free(z_buffer_d);
    free(z_buffer_cor);
    
    return 0;
}

void draw() {
    glBegin(GL_POINTS);
    glColor3f(1, 1, 1);

    int i, j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            if (z_buffer_d[i][j] < INFINITY) {
                glColor3f(z_buffer_cor[i][j][0], z_buffer_cor[i][j][1], z_buffer_cor[i][j][2]);
                glVertex2i(i, j);
            }   
        }
    }
    for (i = 0; i < num_pontos; i++) {
        glVertex2i(pontos_projetados[i][1], pontos_projetados[i][0]);
    }
    glEnd();
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
    pontos =  (float**) calloc(num_pontos, sizeof(float*));
    for(i = 0; i < num_pontos; i++)
    {
        pontos[i] = (float*) calloc(3, sizeof(float));
        fscanf(fp," %f %f %f", &pontos[i][0], &pontos[i][1], &pontos[i][2]);
    }

    triangulos = (int**) calloc(num_triangulos, sizeof(int*));
    for(i = 0; i < num_triangulos; i++)
    {
        triangulos[i] = (int*) calloc(3, sizeof(int));
        fscanf(fp," %d %d %d", &triangulos[i][0], &triangulos[i][1], &triangulos[i][2]);
    }

    fclose (fp);
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
    float aux[3];
    normalizar(N,N);
    proj_vetores(V, N, aux);
    sub_vet(V, aux, V);
    normalizar(V,V);
    //Encontrar U
    prod_vetorial(N, V, U);
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

void iluminar(float V[], float N[], float L[], float ret[]) {
    float difusa[3];
    float ambiental[3];
    float specular[3];
    ambiental[0] = Ia[0]*ka;
    ambiental[1] = Ia[1]*ka;
    ambiental[2] = Ia[2]*ka;
    if (prod_interno(V, N) < 0) {
        mul_escalar(N, -1, N);
    }
    float NL = prod_interno(N, L);
    if(NL>=0){
        difusa[0] = Od[0]*kd*NL*Il[0];
        difusa[1] = Od[1]*kd*NL*Il[1];
        difusa[2] = Od[2]*kd*NL*Il[2];
        float R[3];
        R[0] = 2*NL*N[0];
        R[1] = 2*NL*N[1];
        R[2] = 2*NL*N[2];
        float RV = prod_interno(R, V);
        if(RV>0){
            specular[0] = Il[0]*ks*pow(RV, n);
            specular[1] = Il[1]*ks*pow(RV, n);
            specular[2] = Il[2]*ks*pow(RV, n);
        }
    }

    ret[0] = (ambiental[0]+difusa[0]+specular[0])/255;
    ret[1] = (ambiental[1]+difusa[1]+specular[1])/255;
    ret[2] = (ambiental[2]+difusa[2]+specular[2])/255;
}

void normalizar_triangulos()
{
    int i;
    float aux1[3];
    float aux2[3];
    float aux3[3];
    
    normais_triangulos = (float**) calloc(num_triangulos, sizeof(float*));
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
    float aux[3] = {0,0,0};

    normais_vertices = (float**) calloc(num_pontos, sizeof(float*));
    for(i = 0; i < num_pontos; i++)
    {
        for(j = 0; j < num_triangulos; j++)
        {
            for(k = 0; k < 3; k++)
            {
                if(triangulos[j][k] == (i + 1))
                {
                    sum_vet(aux, normais_triangulos[j], aux);
                }
            }
        }

        normalizar(aux, aux);
        normais_vertices[i] = (float*) calloc(3, sizeof(float));
        normais_vertices[i][0] = aux[0];
        normais_vertices[i][1] = aux[1];
        normais_vertices[i][2] = aux[2];
        aux[0] = 0;
        aux[1] = 0;
        aux[2] = 0;
    }
}

void coord_mundo_para_scc()
{
    int i;
    float** aux = (float**) calloc(num_pontos, sizeof(float*));
    for(i = 0; i < num_pontos; i++)
    {
        aux[i] = (float*) calloc(3, sizeof(float));
        aux[i][0] = U[0]*(pontos[i][0] - C[0]) + U[1]*(pontos[i][1] - C[1]) + U[2]*(pontos[i][2] - C[2]);
        aux[i][1] = V[0]*(pontos[i][0] - C[0]) + V[1]*(pontos[i][1] - C[1]) + V[2]*(pontos[i][2] - C[2]);
        aux[i][2] = N[0]*(pontos[i][0] - C[0]) + N[1]*(pontos[i][1] - C[1]) + N[2]*(pontos[i][2] - C[2]);
        pontos[i][0] = aux[i][0];
        pontos[i][1] = aux[i][1];
        pontos[i][2] = aux[i][2];
    }
    
    free(aux); 
}

void projetar_pontos() {
    pontos_projetados = (int**) calloc(num_pontos, sizeof(int*));
    int i;
    for (i = 0; i < num_pontos; i++) {
        pontos_projetados[i] = (int*) calloc(2, sizeof(int));
        //Esse valor precisa apenas ser multiplicado pela [largura/altura] da tela em que será apresentado e arredondado
        pontos_projetados[i][0] = (int)(((d/hx * pontos[i][0]/pontos[i][2]) + 1)*width/2);
        pontos_projetados[i][1] = (int)((1 - (d/hy * pontos[i][1]/pontos[i][2]))*height/2);
    }
}

void scanline(int** projecao, int** ret) {
    printf("\t\tOrganizando os pontos...");
    int* top = projecao[2];
    int* middle = projecao[1];
    int* bottom = projecao[0];
    printf("OK\n\t\tRealizando scanline do primeiro pedaço...");

    // Calcular o primeiro pedaço dos xmin e xmax
    float a, b;
    int tb_inline = 0, tm_inline = 0;
    if (bottom[0]-top[0] != 0) {
        a = ((float)bottom[1]-top[1])/(bottom[0]-top[0]);
        printf("\t\ta = %f\n", a);
    } else {
        tb_inline = 1;
    }
    if (middle[0]-top[0] != 0) {
        b = ((float)middle[1]-top[1])/(middle[0]-top[0]);
        printf("\t\tb = %f\n", b);
    } else {
        tm_inline = 1;
    }
    int i, j;

    for (i = floor(top[1]), j = 0; i < floor(middle[1]); i++, j++) {
        int x1, x2;
        if (!tb_inline) {
            x1 = floor(top[0] + ((float)j)/a);
        } else {
            x1 = top[0];
        }
        
        if (!tm_inline) {
            x2 = floor(top[0] + ((float)j)/b);
        } else {
            x2 = top[0];
        }

        if (x1 > x2) {
            int tmp = x1;
            x1 = x2;
            x2 = tmp;
        }

        ret[j][0] = x1;
        ret[j][1] = x2;
    }

    printf("OK\n\t\tRealizando scanline do segundo pedaço...");

    // Calcular o segundo pedaço
    int mb_inline = 0;
    if (bottom[0]-middle[0] != 0) {
        b = ((float)bottom[1]-middle[1])/(bottom[0]-middle[0]);
        printf("\t\tb = %f\n", b);
    } else {
        mb_inline = 1;
    }
    for (; i < bottom[1]; i++, j++) {
        int x1, x2;
        if (!tb_inline) {
            x1 = floor(top[0] + ((float)j)/a);
        } else {
            x1 = top[0];
        }
        if (!mb_inline) {
            x2 = floor(middle[0] + ((float)(i-middle[1]))/b);
        } else {
            x2 = middle[0];
        }

        if (x1 > x2) {
            int tmp = x1;
            x1 = x2;
            x2 = tmp;
        }

        ret[j][0] = x1;
        ret[j][1] = x2;
    }

    printf("OK\n");
}

/* Funções do z-buffer */

void init_z_buffer() {
    z_buffer_d = (float**)malloc(height*sizeof(float*));
    z_buffer_cor = (float***)malloc(height*sizeof(float**));

    int i;
    for (i = 0; i < height; i++) {
        z_buffer_d[i] = (float*)malloc(width*sizeof(float));
        z_buffer_cor[i] = (float**)malloc(width*sizeof(float*));

        int j;
        for (j = 0; j < width; j++) {
            z_buffer_d[i][j] = INFINITY;
            z_buffer_cor[i][j] = (float*)malloc(3*sizeof(float));
        }
    }
}

void preencher_z_buffer() {
    int i,j;
    int* projecao[3];
    for (i = 0; i < num_triangulos; i++) {
        printf("\tPegando os pontos do triângulo...");
        // Ordenar os pontos pela coordenada y
        float max = -INFINITY;
        int max_i = 0;
        float min = INFINITY;
        int min_i = 0;
        int* triangulo = triangulos[i];
        for (j = 0; j < 3; j++) {
            if (pontos_projetados[triangulo[j]-1][1] > max) {
                max_i = j;
                max = pontos_projetados[triangulo[j]-1][1];
            } 
       
            if (pontos_projetados[triangulo[j]-1][1] < min) {
                min_i = j;
                min = pontos_projetados[triangulo[j]-1][1];
            }
        }
        int* top = pontos_projetados[triangulo[max_i]-1];
        int* bottom = pontos_projetados[triangulo[min_i]-1];
        int* middle = pontos_projetados[triangulo[0 + 1 + 2 - max_i - min_i]-1];

        projecao[0] = top;
        projecao[1] = middle;
        projecao[2] = bottom;

        printf("OK\n\tRealizando o scanline...\n");
        //printf("\t\tTriangulo:\n\t\t\t(%i, %i)\n\t\t\t(%i, %i)\n\t\t\t(%i, %i)\n", top[0], top[1], middle[0], middle[1], bottom[0], bottom[1]);
        printf("%i %i\n", top[0], top[1]);
        //sleep(5);

        int n_linhas = floor(max) - floor(min);
        int **xminmax;
        printf("\t\tScanline em um triangulo com %i linhas\n", n_linhas);
        if (n_linhas > 0) {
            xminmax = (int**)malloc(n_linhas * sizeof(int*));
            for (j = 0; j < n_linhas; j++) {
                xminmax[j] = (int*)malloc(2*sizeof(int));
            }

            scanline(projecao, xminmax);
        } else {
            continue;
            printf("\t\tCorrigindo e pulando o scanline\n");
            xminmax = (int**)malloc(sizeof(int*));
            xminmax[0] = (int*)malloc(2*sizeof(int));
            xminmax[0][0] = projecao[0][0];
            xminmax[0][1] = projecao[0][0];

            for (j = 0; j < 3; j++) {
                if (projecao[j][0] < xminmax[0][0]) {
                    xminmax[0][0] = projecao[j][0];
                }

                if (projecao[j][0] > xminmax[0][1]) {
                    xminmax[0][1] = projecao[j][0];
                }
            }

            n_linhas = 1;
        }
        
        printf("\t\tResultado do Scanline:\n");
        for (j = 0; j < n_linhas; j++) {
            printf("\t\t\t[%i, %i]\n", xminmax[j][0], xminmax[j][1]);
        }

        printf("                        OK\n\tAcessando o z-buffer...\n");

        int k;
        for (j = 0; j < n_linhas; j++) {
            if (j + bottom[1] < 0 || j + bottom[1] > height) {
                continue;
            }
            for (k = xminmax[j][0]; k <= xminmax[j][1]; k++) {
                if (k < 0 || k > width) {
                    continue;
                }
                printf("\t\tAcessando posição (%i, %i)\n", k, j+bottom[1]);
                printf("\t\tEncontrando coordenadas baricêntricas...");
                int* ponto = (int*)malloc(2*sizeof(int));
                ponto[0] = k; ponto[1] = bottom[1] + j;
                float* coordenadas = (float*)malloc(3*sizeof(float));
                coordenadas_baricentricas(ponto, projecao, coordenadas);
                printf("OK [%f, %f, %f]\n\t\tAproximando o ponto no triangulo...", coordenadas[0], coordenadas[1], coordenadas[2]);
                fflush(stdout);

                float P[3];
                float aux1[3];
                float aux2[3];
                float aux3[3];
                mul_escalar(pontos[triangulo[min_i]-1], coordenadas[0], aux1);
                mul_escalar(pontos[triangulo[max_i]-1], coordenadas[2], aux2);
                sum_vet(aux1, aux2, aux3);
                mul_escalar(pontos[triangulo[3-min_i-max_i]-1], coordenadas[1], aux1);
                sum_vet(aux3, aux1, P);
                printf("OK\n\t\tAtualizando o z-buffer...");

                if (z_buffer_d[bottom[1]+j][k] > P[2]) {
                    // O ponto calculado está mais próximo do que o que está registrado no z-buffer
                    z_buffer_d[bottom[1]+j][k] = P[2];

                    // Mudar a cor salva
                    /// 1. Calcular V, L, N
                    float V[3];
                    float L[3];
                    float N[3];
                    float aux1[3];
                    float aux2[3];
                    float aux3[3];
        
                    mul_escalar(normais_vertices[triangulo[min_i]-1], coordenadas[0], aux1);
                    mul_escalar(normais_vertices[triangulo[max_i]-1], coordenadas[2], aux2);
                    sum_vet(aux1, aux2, aux3);
                    mul_escalar(normais_vertices[triangulo[3-min_i-max_i]-1], coordenadas[1], aux1);
                    sum_vet(aux3, aux1, N);
                    mul_escalar(P, -1, V);
                    sub_vet(Pl, P, L);

                    normalizar(V, V);
                    normalizar(L, L);
                    normalizar(N, N);

                    /// 2. Adicionar a cor ao z-buffer
                    iluminar(V, N, L, z_buffer_cor[bottom[1]+j][k]);
                    printf("\t\tCor do ponto: (R: %f, G: %f, B: %f)\n", z_buffer_cor[bottom[1]+j][k][0], z_buffer_cor[bottom[1]+j][k][1], z_buffer_cor[bottom[1]+j][k][2]);
                }
                printf("OK\n");
            }
        }

        printf("OK\n");
    }

    printf("z-buffer preenchido.\n");
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

void coordenadas_baricentricas(int* ponto, float** triangulo, float* coordenadas) {
    float** sistema = (float**)malloc(3*sizeof(float*));
    
    int i;
    for (i = 0; i < 2; i++) {
        sistema[i] = (float*)malloc(4*sizeof(float));

        int j;
        for (j = 0; j < 3; j++) {
            sistema[i][j] = triangulo[j][i];
        }
        sistema[i][j] = ponto[i];
    }


    sistema[2] = (float*)malloc(4*sizeof(float));
    for (i = 0; i < 4; i++) {
        sistema[2][i] = 1;
    }

    // resolver sistema de equações lineares
    printf("\t\t\tResolvendo o sistema\n");
    resolver_sistema(sistema, 3, 4, coordenadas);
}

void resolver_sistema(float** matriz, int n, int m, float* resultado) {
    printf("\t\t\tEscalonando a matriz do sistema\n");
    escalonar(matriz, n, m);
    
    printf("\t\t\tIsolando as variáveis\n");
    int i;
    for (i = n-1; i >= 0; i--) {
        resultado[i] = matriz[i][m-1];

        int j;
        for (j = i+1; j < n; j++) {
            resultado[i] -= matriz[i][j]*resultado[j];
        }
    }
}

void escalonar(float** matriz, int n, int m) {
    int h = 0; // Inicialização da linha pivô
    int k = 0; // Inicialização da coluna pivô
    
    printf("\t\t\t\tFunção de escalonamento\n");
    while (h < m && k < n) {
        // i_max := argmax(i = h...m, |A(i, k)|)
        int max_i = 0;
        int i;
        for (i = h; i < n; i++) {
            if (fabs(matriz[i][k]) > fabs(matriz[max_i][k])) {
                max_i = i;
            }
        }
        
        if (matriz[max_i][k] == 0) {
            // Não há pivô nessa coluna, passa para a próxima
            k++;
        } else {
            // swap rows(h, i_max)
            float* tmp = matriz[h];
            matriz[h] = matriz[max_i];
            matriz[max_i] = tmp;

            for (i = h+1; i < n; i++) {
                float f = matriz[i][k]/matriz[h][k];

                matriz[i][k] = 0;

                int j;
                for (j = k+1; j < m; j++) {
                    matriz[i][j] = matriz[i][j] - matriz[h][j]*f;
                }
            }

            h++;
            k++;
        }
    }
    printf("\t\t\tEscalonamento completo\n");
}
