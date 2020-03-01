// Copyright (c) 2019. Lloyd T. Elliott

#include "kgen.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdio.h>
#include <time.h>
#include <dirent.h>


#ifndef GIT_VERSION
#  define GIT_VERSION "?"
#endif

#ifndef BUILD_DATE
#  define BUILD_DATE "?"
#endif


#define TREE 1
#define NAIVE 2

static inline long seconds(void) {
  long r = (long)time(NULL);
  return r;
}

int main(int argc, char **argv) {
  opterr = 0;
  int x = -1;
  int method = -1;
  int vflag = 0;
  int dflag = 0;
  int eflag = 0;
  int mflag = 0;
  int flag0 = 0;

  if (argc == 1 || (argc == 2 && (strcmp(argv[1], "--version") == 0))) {
    printf("\n");
    printf("   ┌───────┐   kgen v1 build " GIT_VERSION " " BUILD_DATE ".\n");       
    printf("   │ KGEN  │   Copyright 2019 (c). Lloyd T. Elliott.\n");
    printf("   └───────┘\n");
    printf("\n");

    return 0;
  }

  while ((x = getopt(argc, argv, "0t:vdem:")) != -1) {
    switch (x) {
      case 'v': {
        vflag = 1;
        break;
      }

      case 'e': {
        DOUBLE e = machine_epsilon();
        printf("%e\n", e);
        return 0;
      }

      case 'd': {
        dflag = 1;
        break;
      }

      case '0': {
        flag0 = 1;
        break;
      }
        
      case 'm': {
        mflag = 1;
        if (strcmp(optarg, "tree") == 0) {
          method = TREE;
        } else if (strcmp(optarg, "naive") == 0) {
          method = NAIVE;
        } else {
          error("Unknown method");
        }
        break;
      }

      case '?': {
        error("Argument error");
      }

      default: {
        error("Argument error");
      }
    }
  }

  if (dflag && mflag) {
    error("Argument error");
  }
  if (flag0 && mflag) {
    error("Argument error");
  }
  if (mflag == 0) {
    if (vflag == 1) {
      printf("Using the kgen algorithm.\n\n");
    }
    mflag = 1;
    method = TREE;
  }

  if (dflag && eflag) {
    error("Argument error");
  }
  if (eflag && mflag) {
    error("Argument error");
  }
  if (flag0 && dflag) {
    error("Argument error");
  }
  if (flag0 && eflag) {
    error("Argument error");
  }
  char *first = NULL;
  char *second = NULL;
  char *third = NULL;
  if (flag0 == 1) {
    long t0 = seconds();
    srand((unsigned)time(NULL));
    if (optind > argc - 1) {
      error("More arguments must be provided");
    }
    if (optind < argc - 1) {
      error("Fewer arguments should be provided");
    }
    first = argv[optind];
    INT P1 = find_size(first);
    INT P = sqrt(P1);
    if (P1 != P * P) {
      error("Kinship matrix must be square");
    }

    t_matrix A = load_bin2(first, P, P, P, P);
    t_matrix X = create(P, P, P, P);
    t_matrix Y = clone(A);

    INT d = rand() % P;
    t_matrix V = zeros(P, 1);
    for (int j = 0; j < P; j++) {
      put(V, get(Y, j, d), j, 0);
    }

    long t1 = seconds();
    chol(A);
    long t2 = seconds();
    ddate(&A, d);
    long t3 = seconds();
    udate(&A, d, V);
    long t4 = seconds();
    square(A, 0, &X);
    destroy(V);
    destroy(A);
    DOUBLE f = diff(X, Y);
    destroy(X);
    destroy(Y);
    long t5 = seconds();
    long c1 = t1 - t0;
    long c2 = t2 - t1;
    long c3 = t3 - t2;
    long c4 = t4 - t3;
    long c5 = t5 - t4;
    if (c1 < 1) {
      c1 = 1;
    }
    if (c2 < 1) {
      c2 = 1;
    }
    if (c3 < 1) {
      c3 = 1;
    }
    if (c4 < 1) {
      c4 = 1;
    }
    if (c5 < 1) {
      c5 = 1;
    }
    printf("i c d u l e\n%ld %ld %ld %ld %ld %.10e\n",
      c1,
      c2,
      c3,
      c4,
      c5,
      f);

    exit(0);
  }
  if (dflag == 1) {
    char *dot = NULL;
    if (optind > argc - 2) {
      error("More arguments must be provided");
    }
    if (optind < argc - 2) {
      error("Fewer arguments should be provided");
    }
    first = argv[optind];
    second = argv[optind + 1];
    struct dirent *entry;
    DIR *dp = opendir(second);
    if (dp == NULL) {
      error("Could not open directory");
    }
    int n2 = 0;
    while((entry = readdir(dp))) {
      if (strlen(entry->d_name) <= 2) {
        continue;
      }
      dot = strrchr(entry->d_name, '.');
      if (dot == entry->d_name || dot == &(entry->d_name[1])) {
        continue;
      }
      n2++;
    }
    closedir(dp);
    dp = opendir(first);
    if (dp == NULL) {
      error("Could not open directory");
    }
    int n1 = 0;
    while((entry = readdir(dp))) {
      if (strlen(entry->d_name) <= 2) {
        continue;
      }
      dot = strrchr(entry->d_name, '.');
      if (dot == entry->d_name || dot == &(entry->d_name[1])) {
        continue;
      }
      n1++;
    }

    closedir(dp);
    if (n1 != n2) {
      error("Number of phenotypes differ");
    }
    dp = opendir(first);
    if (dp == NULL) {
      error("Could not open directory");
    }
    while ((entry = readdir(dp))) {
      if (strlen(entry->d_name) <= 2) {
        continue;
      }
      dot = strrchr(entry->d_name, '.');
      if (dot == entry->d_name || dot == &(entry->d_name[1])) {
        continue;
      }
      n1++;
      size_t needed1 = snprintf(NULL, 0, "%s/%s", first, entry->d_name);
      char *fname1 = NULL;
      if ((fname1 = (char *)malloc(sizeof(char) * (1 + needed1))) == NULL) {
        error("Memory exhausted");
      }
      sprintf(fname1, "%s/%s", first, entry->d_name);
      INT P = find_size(fname1);
      t_matrix A = load_bin3(fname1, P, 1, P, 1);
      free(fname1);
      size_t needed2 = snprintf(NULL, 0, "%s/%s", second, entry->d_name);
      char *fname2 = NULL;
      if ((fname2 = (char *)malloc(sizeof(char) * (1 + needed2))) == NULL) {
        error("Memory exhausted");
      }
      sprintf(fname2, "%s/%s", second, entry->d_name);
      t_matrix B = load_bin3(fname2, P, 1, P, 1);
      free(fname2);
      DOUBLE f = diff(A, B);
      printf("%.10e\n", f);
      destroy(A);
      destroy(B);
    }

    closedir(dp);
    return 0;
  }
  
  if (optind > argc - 3) {
    error("More arguments must be provided");
  }
  
  if (optind + 3 < argc) {
    error("Fewer arguments should be provided");
  }

  first = argv[optind];
  second = argv[optind + 1];
  third = argv[optind + 2];
  struct stat sb;
  int e = stat(third, &sb);
  if (e == 0) {
    if (!S_ISDIR(sb.st_mode)) {
      error("Cannot create output directory");
    }
  } else {
    if (errno = ENOENT) {
      e = mkdir(third, S_IRWXU);
      if (e != 0) {
       error("Could not create output directory");
      }
    } else {
      error("i/o error");
    }
  }

  switch (method) {
    case NAIVE: {
      long ta = seconds();
      char *dot = NULL;
      INT N = 0;
      INT D = 0;
      dot = strrchr(second, '.');
      if (dot && !strcmp(dot, ".bin")) {
        INT P1 = find_size(second);
        INT P = sqrt(P1);
        if (P1 != P * P) {
          error("Kinship matrix must be square");
        }
        N = P;
      }

      t_matrix X;
      dot = strrchr(first, '.');
      if ((dot != NULL) && !strcmp(dot, ".bin")) {
        INT size = find_size(first);
        if (N == 0) {
          error("Binary files must be provided");
        }
        D = size / N;
        if (size != N * D) {
          error("Phenotype matrix must have N rows");
        }
        X = load_bin2(first, N, D, N, D);
      } else {
        X = load(first, 1);
        if (N > 0 && N != X.N) {
          error("Dimension of phenotype and kinship mismatch");
        }
        N = X.N;
        D = X.M;
      }

      is_nan(X);
      if (vflag) {
        printf("Loading %lld phenotypes and %lld subjects.\n", D, N);
        printf("Number of missing entries:");
        int nn = 0;
        for (INT i = 0; i < D; i++) {
          for (INT k = 0; k < N; k++) {
            if (get(X, k, i) > 0.5) {
              nn += 1;
            }
          }
        }
        printf(" %d\n", nn);
      }
      t_matrix A = zeros(D, D);
      for (INT i = 1; i < D; i++) {
        for (INT j = 0; j < i; j++) {
           put(A, 1, i, j);
        }
      }

      for (INT i = 1; i < D; i++) {
        for (INT j = 0; j < i; j++) {
          for (INT k = 0; k < N; k++) {
            if ((get(X, k, i) > 0.5 && get(X, k, j) < 0.5) ||
              (get(X, k, i) < 0.5 && get(X, k, j) > 0.5)) {

              put(A, get(A, i, j) + 1, i, j);
            } 
          }
        }
      }

      for (INT i = 1; i < D; i++) {
        for (INT j = 0; j < i; j++) {
          put(A, get(A, i, j), j, i);
        }
      }

      t_graph G = create_graph(A);
      destroy(A);
      dot = strrchr(second, '.');
      if ((dot != NULL) && !strcmp(dot, ".bin")) {
        A = load_bin2(second, N, N, N, N);
      } else {
        A = load(second, 0);
        if (A.M != N || A.N != N) {
          error("Mismatch between dimension of phenotype and kinship");
        }
      }
      
      INT *which = NULL;
      if ((which = (INT *)malloc((N + 1) * sizeof(INT))) == NULL) {
        error("Memory exhausted");
      }
      for (int i = 0; i < D; i++) {
        put(A, get(A, i, i) + 1e-3, i, i);
      }
      for (int id = 0; id < D; id++) {
        if (vflag) {
          fflush(stdout);
        }
        int found = 0;
        if (id > 0) {
          for (int id0 = 0; id0 < id; id0++) {
            if (get_weight(G, id0, id) < 1.5) {
              found = 1;
              int n1 = id0;
              int n2 = id;
              int needed1 = snprintf(NULL, 0, "%s/L%05d.bin", third, n1);
              char *fname1 = (char *)malloc(sizeof(char) * (1 + needed1));
              if (fname1 == NULL) {
                error("Memory exhausted");
              }
              sprintf(fname1, "%s/L%05d.bin", third, n1);
              int needed2 = snprintf(NULL, 0, "%s/L%05d.bin", third, n2);
              char *fname2 = (char *)malloc(sizeof(char) * (1 + needed2));
              if (fname2 == NULL) {
                error("Memory exhausted");
              }
              sprintf(fname2, "%s/L%05d.bin", third, n2);
              if (vflag) {
                printf("Copying decomposition %d - %d\n", id0, id);
              }
              cp(fname1, fname2);
              free(fname1);
              free(fname2);
              break;
            }
          }

          if (found) {
            continue;
          }
        }

        int k = 0;
        for (INT i = 0; i < N; i++) {
          if (get(X, i, id) == 1) {
            which[k] = i;
            k++;
          }
        }

        which[k] = -1;
        t_matrix K = clone(A);
        remove_square(&K, k, which);
        if (vflag) {
          printf("Performing decomposition %d\n", id);
        }
        chol(K);
        size_t needed = snprintf(NULL, 0, "%s/L%05d.bin", third, id);
        char *fname = NULL;
        if ((fname = (char *)malloc(sizeof(char) * (1 + needed))) == NULL) {
          error("Memory exhausted");
        }
        sprintf(fname, "%s/L%05d.bin", third, id);
        save_bin3(K, fname);
        free(fname);
        destroy(K);
      }

      destroy(A);
      destroy(X);
      destroy(G);
      free(which);
      long tb = seconds();
      long total = ta - tb;
      if (total < 1) {
        total = 1;
      }
      printf("t\n%ld\n", total);
      break;
    }

    case TREE: {
      long ta = seconds();
      size_t needed = 0;
      char *fname = NULL;
      char *dot = NULL;
      INT N = 0;
      INT D = 0;
      dot = strrchr(second, '.');
      if (dot && !strcmp(dot, ".bin")) {
        INT P1 = find_size(second);
        INT P = sqrt(P1);
        if (P1 != P * P) {
          error("Kinship matrix must be square");
        }
        N = P;
      } 

      t_matrix A;
      if (first == NULL) {
        error("Phenotype file must be provided");
      }

      t_matrix X;
      dot = strrchr(first, '.');
      if ((dot != NULL) && !strcmp(dot, ".bin")) {
        INT size = find_size(first);
        if (N == 0) {
          error("Binary files must be provided");
        }

        D = size / N;
        if (size != N * D) {
          error("Phenotype matrix must have N rows");
        }

        X = load_bin2(first, N, D, N, D);
      } else {
        X = load(first, 1);
        if (N > 0 && N != X.N) {
          error("Dimension of phenotype and kinship mismatch");
        }
        N = X.N;
        D = X.M;

      }

      if (N == 0) {
        error("Number of subjects unknown");
      }

      if (D == 0) {
        error("Number of dimensions unknown");
      }

      is_nan(X);
      if (vflag) {
        printf("Loading %lld phenotypes and %lld subjects.\n", D, N);
        printf("Number of missing entries:");
        int nn = 0;
        for (INT i = 0; i < D; i++) {
          for (INT k = 0; k < N; k++) {
            if (get(X, k, i) > 0.5) {
              nn += 1;
            }
          }
        }
        printf(" %d\n", nn);
      }
      A = zeros(D, D);
      for (INT i = 1; i < D; i++) {
        for (INT j = 0; j < i; j++) {
           put(A, 1, i, j);
        }
      }

      for (INT i = 1; i < D; i++) {
        for (INT j = 0; j < i; j++) {
          for (INT k = 0; k < N; k++) {
            if ((get(X, k, i) > 0.5 && get(X, k, j) < 0.5) ||
              (get(X, k, i) < 0.5 && get(X, k, j) > 0.5)) {

              put(A, get(A, i, j) + 1, i, j);
            } 
          }
        }
      }

      for (INT i = 1; i < D; i++) {
        for (INT j = 0; j < i; j++) {
          put(A, get(A, i, j), j, i);
        }
      }

      t_graph G = create_graph(A);
      minspan(G);
      destroy(A);
      if (second == NULL) {
        error("Kinship matrix must be provided");
      }

      dot = strrchr(second, '.');
      if ((dot != NULL) && !strcmp(dot, ".bin")) {
        A = load_bin2(second, N, N, N, N);
      } else {
        A = load(second, 0);
        if (A.M != N || A.N != N) {
          error("Mismatch between dimension of phenotype and kinship");
        }
      }

      for (int i = 0; i < D; i++) {
        put(A, get(A, i, i) + 1e-3, i, i);
      }

      t_matrix K = clone(A);
      INT *which = NULL;
      if ((which = (INT *)malloc((N+1) * sizeof(INT))) == NULL) {
        error("Memory exhausted");
      }

      int **ids = NULL;
      if ((ids = (int **)malloc(D * sizeof(int *))) == NULL) {
        error("Memory exhausted");
      }

      for (int j = 0; j < D; j++) {
        if ((ids[j] = (int *)malloc(N * sizeof(int))) == NULL) {
          error("Memory exhausted");
        }
      }

      int *lengths = NULL;
      if ((lengths = (int *)malloc(D * sizeof(int))) == NULL) {
        error("Memory exhausted");
      }

      bzero(lengths, D * sizeof(int));
      for (int j = 0; j < D; j++) {
        int m = 0;
        for (int i = 0; i < N; i++) {
          ids[j][m] = -1;
        }

        for (int i = 0; i < N; i++) {
          if (iszero(get(X, i, j))) {
            lengths[j]++;
            ids[j][m] = i;
            m++;
          }
        }
      }

      int id = 0;
      int mm = lengths[0];
      int k = 0;
      for (INT j = 0; j < D; j++) {
        if (lengths[j] < mm) {
          id = j;
          mm = lengths[j];
        }
      }

      for (INT i = 0; i < N; i++) {
        if (get(X, i, id) == 1) {
          which[k] = i;
          k++;
        }
      }

      which[k] = -1;
      remove_square(&A, k, which);
      free(which);
      chol(A);
      needed = snprintf(NULL, 0, "%s/L%05d.bin", third, id);
      fname = NULL;
      if ((fname = (char *)malloc(sizeof(char) * (1 + needed))) == NULL) {
        error("Memory exhausted");
      }

      sprintf(fname, "%s/L%05d.bin", third, id);
      save_bin3(A, fname);
      free(fname);
      
      int *done = NULL;
      if ((done = (int *)malloc(D * sizeof(int))) == NULL) {
        error("Memory exhausted");
      }
      bzero(done, D * sizeof(int));
      t_edges edges = get_edges(G);
      for (int i = 0; i < D - 1; i++) { // Find root
        if (edges.e2[i] == id) {
          SWAP(edges.e1[i], edges.e2[i]);
        }
      }

      int first = -1;
      for (int i = 0; i < D - 1; i++) {
        if (edges.e1[i] == id && first < 0) {
          first = i;
          break;
        }
      }

      if (first < 0) {
        error("Could not find root");
      }
      
      SWAP(edges.e1[first], edges.e1[0]);
      SWAP(edges.e2[first], edges.e2[0]);
      done[edges.e1[0]] = 1;
      done[edges.e2[0]] = 1;
      int current = id;
      for (int i = 1; i < D - 1; i++) { // Breadth first
        int found = -1;
        for (int j = i; j < D - 1 && found == -1; j++) {
          if (done[edges.e1[j]] && edges.e1[j] == current) {
            found = j;
          } else if (done[edges.e2[j]] && edges.e2[j] == current) {
            SWAP(edges.e1[j], edges.e2[j]);
            found = j;
          }
        }

        if (found == -1) {
          for (int j = i; j < D - 1 && found == -1; j++) {
            if (done[edges.e1[j]]) {
              found = j;
              current = edges.e1[j];
            } else if (done[edges.e2[j]]) {
              SWAP(edges.e1[j], edges.e2[j]);
              found = j;
              current = edges.e1[j];
            }
          }
        }

        if (found == -1) {
          error("Could not find edge");
        }
        if (found > 0) {
          done[edges.e1[found]] = 1;
          done[edges.e2[found]] = 1;
          SWAP(edges.e1[i], edges.e1[found]);
          SWAP(edges.e2[i], edges.e2[found]);
          assert(current == edges.e1[i]);
        }
      }

      if (vflag) {
        printf("Considering path:\n");
        for (int i = 0; i < D - 1; i++) {
          printf("  %d - %d (%d)\n",
            edges.e1[i],
            edges.e2[i],
            (int)get_weight(G, edges.e1[i], edges.e2[i]) - 1);

        }

        printf("\nTotal path weight:");
        double w = 0.0;
        for (int i = 0; i < D - 1; i++) {
          w += get_weight(G, edges.e1[i], edges.e2[i]) - 1;
        }
        printf(" %f\n", w);
      }

      free(done);
      assert(edges.e1[0] == id);
      int prev = id;
      for (int i = 0; i < D - 1; i++) {
        if (vflag) {
          fflush(stdout);
        }
        t_matrix L2;
        int n1 = edges.e1[i];
        int n2 = edges.e2[i];
        
        if (get_weight(G, n1, n2) < 1.5) {
          int needed1 = snprintf(NULL, 0, "%s/L%05d.bin", third, n1);
          char *fname1 = NULL;
          if ((fname1 = (char *)malloc(sizeof(char) * (1 + needed1))) == NULL) {
            error("Memory exhausted");
          }
          sprintf(fname1, "%s/L%05d.bin", third, n1);
          int needed2 = snprintf(NULL, 0, "%s/L%05d.bin", third, n2);
          char *fname2 = NULL;
          if ((fname2 = (char *)malloc(sizeof(char) * (1 + needed2))) == NULL) {
            error("Memory exhausted");
          }
          sprintf(fname2, "%s/L%05d.bin", third, n2);
          cp(fname1, fname2);
          free(fname1);
          free(fname2);
          if (vflag) {
            printf("Copying edge %d - %d.\n", n1, n2);
          }

          if (prev != edges.e1[i]) {
            needed = snprintf(NULL, 0, "%s/L%05d.bin", third, n1);
            fname = NULL;
            if ((fname = (char *)malloc(sizeof(char) * (1 + needed))) == NULL) {
              error("Memory exhausted");
            }
            sprintf(fname, "%s/L%05d.bin", third, n1);
            destroy(A);
            A = load_bin3(fname, lengths[n1], lengths[n1], N, N);
            free(fname);
            prev = edges.e1[i];
            if (vflag) {
              printf(" Reading decomposition %d from disk\n", prev);
            }
          }
          continue;
        } else if (get_weight(G, n1, n2) > 100) {
          if (vflag) {
            printf(" Computing phenotype %d.\n", n2);
          }
          t_matrix L2 = clone(K);
          INT *which = NULL;
          if ((which = (INT *)malloc((N+1) * sizeof(INT))) == NULL) {
            error("Memory exhausted");
          }

          int k = 0;
          for (INT i = 0; i < N; i++) {
            if (get(X, i, n2) == 1) {
              which[k] = i;
              k++;
            }
          }

          which[k] = -1;
          remove_square(&L2, k, which);
          free(which);
          chol(L2);
          needed = snprintf(NULL, 0, "%s/L%05d.bin", third, n2);
          fname = NULL;
          if ((fname = (char *)malloc(sizeof(char) * (1 + needed))) == NULL) {
            error("Memory exhausted");
          }
          sprintf(fname, "%s/L%05d.bin", third, n2);
          save_bin3(L2, fname);
          free(fname);
          destroy(L2);
          if (prev != edges.e1[i]) {
            needed = snprintf(NULL, 0, "%s/L%05d.bin", third, n1);
            fname = NULL;
            if ((fname = (char *)malloc(sizeof(char) * (1 + needed))) == NULL) {
              error("Memory exhausted");
            }
            sprintf(fname, "%s/L%05d.bin", third, n1);
            destroy(A);
            A = load_bin3(fname, lengths[n1], lengths[n1], N, N);
            free(fname);
            prev = edges.e1[i];
            if (vflag) {
              printf(" Reading decomposition %d from disk\n", prev);
            }
          }
          printf("F\n");
          continue;
        }

        if (vflag) {
          printf("Considering edge %d - %d (%d)\n",
            n1,
            n2,
            (int)get_weight(G, n1, n2) - 1);

        }

        if (prev == edges.e1[i]) {
          L2 = clone(A);
        } else {
          needed = snprintf(NULL, 0, "%s/L%05d.bin", third, n1);
          fname = NULL;
          if ((fname = (char *)malloc(sizeof(char) * (1 + needed))) == NULL) {
            error("Memory exhausted");
          }
          sprintf(fname, "%s/L%05d.bin", third, n1);
          L2 = load_bin3(fname, lengths[n1], lengths[n1], N, N);
          free(fname);
          destroy(A);
          A = clone(L2);
          prev = edges.e1[i];
          if (vflag) {
            printf(" Reading decomposition %d from disk\n", prev);
          }
        }

        if (n1 == n2) {
          error("Edges are wrong");
        }
        int *zs = NULL;
        if ((zs = (int *)malloc(N * sizeof(int))) == NULL) {
          error("Memory exhausted");
        }
        memcpy(zs, ids[n1], lengths[n1] * sizeof(int));
        int length = lengths[n1];
        int z1 = lengths[n1] - 1;
        int z0 = z1;
        int z2 = lengths[n2] - 1;

        for (int n = N - 1; n >= 0; n--) {
          if (get(X, n, n1) == 0 && get(X, n, n2) == 0) {
            z0 -= 1;
            z1 -= 1;
            z2 -= 1;
          } else if (get(X, n, n1) == 1 && get(X, n, n2) == 0) {
            if (vflag) {
              printf("  Inserting subject %d\n", n);
            }
            z2 -= 1;
            if (length - z0 + 1 > 0) {
              memmove(&(zs[z0 + 2]),
                &(zs[z0 + 1]),
                sizeof(int) * (length - (z0 + 1)));

            }

            zs[z0 + 1] = n;
            length += 1;
            t_matrix V = zeros(length, 1);
            for (int j = 0; j < length; j++) {
              put(V, get(K, zs[j], n), j, 0);
            }
            udate(&L2, z0 + 1, V);
            destroy(V);
          } else if (get(X, n, n1) == 0 && get(X, n, n2) == 1) {
            if (vflag) {
              printf("  Deleting subject %d\n", n);
            }

            z1 -= 1;
            if (length - z0 + 1 > 0) {
              memmove(&(zs[z0]),
                &(zs[z0 + 1]),
                sizeof(int) * (length - (z0 + 1)));

            }

            length -= 1;
            if (length == 0) {
              error("Refusing to operate as all samples must be removed");
            }

            ddate(&L2, z0);
            z0 -= 1;

          } else if (get(X, n, n1) == 1 && get(X, n, n2) == 1) {

          }
        }

        int *zs2 = (int *)malloc(sizeof(int) * length);
        int kk = 0;
        for (int i = 0; i < N; i++) {
          if (get(X, i, n2) == 0) {
            zs2[kk] = i;
            kk++;
          }

          if (kk > length) {
            error("Assertion failed");
          }
        } 

        for (int i = 0; i < length; i++) {
          if (zs2[i] != zs[i]) {
            error("Assertion failed");
          }
        }

        free(zs2);
        free(zs);
        needed = snprintf(NULL, 0, "%s/L%05d.bin", third, n2);
        fname = NULL;
        if ((fname = (char *)malloc(sizeof(char) * (1 + needed))) == NULL) {
          error("Memory exhausted");
        }
        sprintf(fname, "%s/L%05d.bin", third, n2);
        save_bin3(L2, fname);
        free(fname);
        destroy(L2);
      }

      if (vflag) {
        printf("\n");
      }
      for (int j = 0; j < D; j++) {
        free(ids[j]);
      }
      free(ids);
      free(lengths);
      destroy(G);
      destroy(A);
      destroy(K);
      destroy(X);
      destroy(edges);
      long tb = seconds();
      long total = ta - tb;
      if (total < 1) {
        total = 1;
      }
      printf("t\n%ld\n", total);
      break;
    }

    default: {
      error("Unknown method");
    }
  }

  exit(0);
}
