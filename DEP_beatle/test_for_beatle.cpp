// This code represents model predictive controller functionality for driving envelope protection in C++
// It utilizes qpOASES solver
// Author: Denis Efremov from CTU in Prague

#include <qpOASES.hpp>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

USING_NAMESPACE_QPOASES

void makeSparseFromDense(real_t *denseMatrix,
         int_t colNumber,
         int_t rowNumber,
         int_t nnz,
         int_t *rowidx,
         int_t *colidx,
         real_t *d) {
    int_t pos;
    real_t data;
    rowidx[0] = 0;
    colidx[0] = 0;
    pos = 0;

    for (int_t col{0}; col < colNumber; col++) {
        for (int_t row{0}; row < rowNumber; row++) {
            data = denseMatrix[row*colNumber+col]; // trying change, because I'm obtaining transposed matrix
            if (data != 0.0) {
                rowidx[pos] = row;
                d[pos] = data;
                pos++;
            }
        }
        colidx[col + 1] = pos;
    }
}

// autogenerated function from matlab
// it defines matrices for the OCP before each call of the solver
void defineOPC(const real_t input[9], const real_t params[9],
               real_t Q[1089], real_t c[33],
               real_t Asolver[3069], real_t lb[33],
               real_t ub[33], real_t lbA[93], real_t ubA[93]){
    static const double y[16]{40.0, 0.0, 0.0,  0.0, 0.0, 0.01, 0.0, 0.0,
                              0.0,  0.0, 40.0, 0.0, 0.0, 0.0,  0.0, 0.01};
    static const double dv2[6]{200.0, 0.6, 200.0, 0.6, 200.0, 0.6};
    static const double dv[5]{-1.3333333333333335, 3.3333333333333335,
                              -1.3333333333333335, 3.3333333333333335, 0.0};
    static const double dv1[5]{-1.3333333333333335, 3.3333333333333335,
                               -1.3333333333333335, 3.3333333333333335, 0.0};
    static const double b_y[2]{100.0, 0.3};
    double A[2871];
    double Aeq[198];
    double b[87];
    double Q_tmp[36];
    double b_d[36];
    double d[36];
    double b_statesFL[10];
    double inputsFL_tmp[8];
    double statesFL[8];
    double beq[6];
    double Ad[4];
    double Bd[4];
    double slew[2];
    double b_params_idx_0_tmp;
    double b_params_idx_1_tmp;
    double lr;
    double omega_slew;
    double params_idx_0_tmp;
    double params_idx_1_tmp;
    double params_idx_1_tmp_tmp;
    double v;
    double wf_now;
    double x1_idx_0;
    double x1_idx_1;
    int A_tmp;
    int A_tmp_tmp;
    int b_A_tmp;
    int b_j;
    int i;
    int j;
    int j_tmp;
    signed char Ad_tmp[4];
    //  vehicle parameters
    //  GC-front distance
    lr = params[1];
    //  GC-rear distance
    //  mass
    //  inertia
    //  half of the axle
    //  wheel radius
    //  cornering stiffness front
    //  cornering stiffness rear
    //  input parsing
    x1_idx_0 = input[0];
    x1_idx_1 = input[1];
    v = input[2];
    wf_now = input[3];
    //      %% define the system
    Ad_tmp[1] = 0;
    Ad_tmp[2] = 0;
    Ad_tmp[0] = 1;
    Ad_tmp[3] = 1;
    params_idx_0_tmp =
            params[8] * (-(params[6] + params[7]) / params[2] / input[2]);
    params_idx_1_tmp_tmp = params[0] * params[6];
    params_idx_1_tmp = params[1] * params[7] - params_idx_1_tmp_tmp;
    b_params_idx_1_tmp =
            params[8] * (params_idx_1_tmp / params[2] / input[2] / input[2] - 1.0);
    b_params_idx_0_tmp = params[8] * (params_idx_1_tmp / params[3]);
    params_idx_1_tmp = params[8] * (-(params[0] * params[0] * params[6] +
                                      params[7] * (params[1] * params[1])) /
                                    params[3] / input[2]);
    Ad[0] = params_idx_0_tmp + 1.0;
    Ad[1] = b_params_idx_0_tmp;
    Ad[2] = b_params_idx_1_tmp;
    Ad[3] = params_idx_1_tmp + 1.0;
    Bd[0] = params[8] * (params[6] / params[2] / input[2]);
    Bd[2] = params[8] * 0.0;
    Bd[1] = params[8] * (params_idx_1_tmp_tmp / params[3]);
    Bd[3] = params[8] * 0.0;
    //  mpc stuff
    //  number of states
    //  number of inputs
    //  Weights
    //  100 20
    // 0.008;
    //  slack weight % 1e4
    // envelope weight
    //  0.08 0,09 0,5
    //      if(thr_pedal > 0.1)
    //          omega_track = 0.5;
    //  %         Qdu(2) = 0.5;
    //      end
    //  MPC data
    //  changable envelope restrictions
    //      v = sdpvar(1,1); % velocity at the moment
    //  slew restrictions
    // rad/s
    omega_slew = 400.0 * params[8];
    //  150 % 600 is the best for braking
    slew[0] = 2.0943951023931953 * params[8];
    //  u11,u12,...,u1N
    //  u21,u22,...,u2N
    //  x12,...,x1N+1
    //  x22,...,x2N+1
    //  sAR2,sAR3,...sARN+1 % in step 1 AR cannot be optimized
    //  sD1,sD2,...,sDN
    //  sW1,sW2,...,sWN
    //  sSFL1,sSFL2,...sSFLN+1
    //  sSFR1,sSFR2,...sSFRN+1
    //  c11,c12,...,c1N
    //  c21,c22,...,c2N
    //  count of abs in objective function
    // (states+inputs)*N+slacks+abs_num*inputs
    std::memset(&Q[0], 0, 1089U * sizeof(double));
    std::memset(&c[0], 0, 33U * sizeof(double));
    std::memset(&A[0], 0, 2871U * sizeof(double));
    //  inequality constraints (min/max, aR and its slacks, du and its slacks,
    //  sigmaFL and sigmaFR and its slacks abs)
    std::memset(&b[0], 0, 87U * sizeof(double));
    std::memset(&Aeq[0], 0, 198U * sizeof(double));
    // equality constraints (dynamics)
    for (i = 0; i < 6; i++) {
        beq[i] = 0.0;
    }
    //  boundaries for variables (seems to be unused in yalmip)
    //     %% min/max inputs
    for (j = 0; j < 3; j++) {
        A_tmp_tmp = j << 2;
        j_tmp = j << 1;
        A_tmp = A_tmp_tmp + 87 * j_tmp;
        A[A_tmp] = -1.0;
        A[A_tmp + 1] = 1.0;
        b[A_tmp_tmp] = 0.65;
        b[A_tmp_tmp + 1] = 0.65;
        A_tmp = A_tmp_tmp + 87 * (j_tmp + 1);
        A[A_tmp + 2] = -1.0;
        A[A_tmp + 3] = 1.0;
        params_idx_1_tmp_tmp = (static_cast<double>(j) + 1.0) * omega_slew;
        b[A_tmp_tmp + 2] = -(wf_now - params_idx_1_tmp_tmp);
        b[A_tmp_tmp + 3] = wf_now + params_idx_1_tmp_tmp;
    }
    //     %% aR envelope and slack constraints
    for (j = 0; j < 3; j++) {
        //  ar
        j_tmp = j * 3 + 13;
        b_j = (j << 1) + 6;
        A_tmp = j_tmp + 87 * b_j;
        A[A_tmp - 1] = -1.0;
        b_A_tmp = j_tmp + 87 * (b_j + 1);
        A[b_A_tmp - 1] = lr / v;
        A[A_tmp] = 1.0;
        A[b_A_tmp] = -lr / v;
        A[A_tmp + 1] = 0.0;
        A[b_A_tmp + 1] = 0.0;
        j_tmp = j * 3;
        A_tmp = j * 3 + 87 * (j + 12);
        A[A_tmp + 12] = -1.0;
        b[j_tmp + 12] = 0.4;
        A[A_tmp + 13] = -1.0;
        b[j_tmp + 13] = 0.4;
        A[A_tmp + 14] = -1.0;
        b[j_tmp + 14] = 0.0;
        Q[(j + 33 * (j + 12)) + 12] = 2.0E+6;
    }
    //     %% slew protection
    for (j = 0; j < 3; j++) {
        if (j + 1 == 1) {
            b[21] = -input[7] + slew[0];
            b[23] = slew[0] + input[7];
            A[21] = -1.0;
            A[1326] = -1.0;
            A[23] = 1.0;
            A[1328] = -1.0;
            A[1330] = -1.0;
            A[22] = 0.0;
            A[1327] = 0.0;
            A[24] = 0.0;
            A[1329] = 0.0;
            A[1331] = 0.0;
            b[22] = -input[8] + omega_slew;
            b[24] = omega_slew + input[8];
            A[108] = 0.0;
            A[1413] = 0.0;
            A[110] = 0.0;
            A[1415] = 0.0;
            A[1417] = 0.0;
            A[109] = -1.0;
            A[1414] = -1.0;
            A[111] = 1.0;
            A[1416] = -1.0;
            A[1418] = -1.0;
        } else {
            int b_A_tmp_tmp;
            int b_j_tmp;
            j_tmp = j << 1;
            b_j_tmp = j_tmp * 3;
            A_tmp = b_j_tmp + 87 * (j_tmp - 2);
            A[A_tmp + 21] = 1.0;
            A[A_tmp + 22] = 0.0;
            b_A_tmp_tmp = 87 * (j_tmp - 1);
            b_A_tmp = b_j_tmp + b_A_tmp_tmp;
            A[b_A_tmp + 21] = 0.0;
            A[b_A_tmp + 22] = 1.0;
            b_j = b_j_tmp + 87 * j_tmp;
            A[b_j + 21] = -1.0;
            A[b_j + 22] = 0.0;
            i = b_j_tmp + 87 * (j_tmp + 1);
            A[i + 21] = 0.0;
            A[i + 22] = -1.0;
            A_tmp_tmp = b_j_tmp + 87 * (j_tmp + 15);
            A[A_tmp_tmp + 21] = -1.0;
            A[A_tmp_tmp + 22] = 0.0;
            b[b_j_tmp + 21] = slew[0];
            j_tmp = b_j_tmp + 87 * (j_tmp + 16);
            A[j_tmp + 21] = 0.0;
            A[j_tmp + 22] = -1.0;
            b[b_j_tmp + 22] = omega_slew;
            A[A_tmp + 23] = -1.0;
            A[A_tmp + 24] = 0.0;
            A_tmp = b_j_tmp + b_A_tmp_tmp;
            A[A_tmp + 23] = 0.0;
            A[A_tmp + 24] = -1.0;
            A[b_j + 23] = 1.0;
            A[b_j + 24] = 0.0;
            A[i + 23] = 0.0;
            A[i + 24] = 1.0;
            A[A_tmp_tmp + 23] = -1.0;
            A[A_tmp_tmp + 24] = 0.0;
            b[b_j_tmp + 23] = slew[0];
            A[j_tmp + 23] = 0.0;
            A[j_tmp + 24] = -1.0;
            b[b_j_tmp + 24] = omega_slew;
            A[A_tmp_tmp + 25] = -1.0;
            A[A_tmp_tmp + 26] = 0.0;
            A[j_tmp + 25] = 0.0;
            A[j_tmp + 26] = -1.0;
        }
        j_tmp = j << 1;
        b_j = j << 1;
        i = j_tmp + 33 * (b_j + 15);
        Q[i + 15] = 2000.0;
        Q[i + 16] = 0.0;
        i = j_tmp + 33 * (b_j + 16);
        Q[i + 15] = 0.0;
        Q[i + 16] = 2000.0;
    }
    wf_now = 39.0;
    //     %% sigmaFL envelope and slack constrains
    inputsFL_tmp[0] = -1.6556556940273772;
    params_idx_1_tmp_tmp = -2.3333333333333335 * params[5] / input[2];
    inputsFL_tmp[4] = params_idx_1_tmp_tmp;
    inputsFL_tmp[1] = 1.6556556940273772;
    omega_slew = 2.3333333333333335 * params[5] / input[2];
    inputsFL_tmp[5] = omega_slew;
    inputsFL_tmp[2] = 1.6556556940273772;
    inputsFL_tmp[6] = params_idx_1_tmp_tmp;
    inputsFL_tmp[3] = -1.6556556940273772;
    inputsFL_tmp[7] = omega_slew;
    statesFL[0] = 1.6556556940273772;
    params_idx_1_tmp_tmp =
            (2.3333333333333335 * params[4] + 1.6556556940273772 * params[0]) /
            input[2];
    statesFL[4] = params_idx_1_tmp_tmp;
    statesFL[1] = -1.6556556940273772;
    omega_slew =
            (-2.3333333333333335 * params[4] - 1.6556556940273772 * params[0]) /
            input[2];
    statesFL[5] = omega_slew;
    statesFL[2] = -1.6556556940273772;
    lr = (2.3333333333333335 * params[4] - 1.6556556940273772 * params[0]) /
         input[2];
    statesFL[6] = lr;
    statesFL[3] = 1.6556556940273772;
    v = (-2.3333333333333335 * params[4] + 1.6556556940273772 * params[0]) /
        input[2];
    statesFL[7] = v;
    i = 21;
    for (b_j = 0; b_j < 3; b_j++) {
        i = b_j + 21;
        if (b_j + 1 == 1) {
            j_tmp = b_j << 1;
            for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
                A_tmp = A_tmp_tmp << 2;
                b_A_tmp = 87 * (A_tmp_tmp + j_tmp);
                A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp];
                A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 1];
                A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 2];
                A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 3];
            }
            A[(static_cast<int>(wf_now + 5.0) + 87 * j_tmp) - 1] = 0.0;
            A[(static_cast<int>(wf_now + 5.0) + 87 * (j_tmp + 1)) - 1] = 0.0;
            for (A_tmp_tmp = 0; A_tmp_tmp < 5; A_tmp_tmp++) {
                A[(static_cast<int>(wf_now + (static_cast<double>(A_tmp_tmp) + 1.0)) +
                   87 * (b_j + 21)) -
                  1] = -1.0;
            }
            for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
                j_tmp = A_tmp_tmp << 2;
                b_statesFL[5 * A_tmp_tmp] = statesFL[j_tmp];
                b_statesFL[5 * A_tmp_tmp + 1] = statesFL[j_tmp + 1];
                b_statesFL[5 * A_tmp_tmp + 2] = statesFL[j_tmp + 2];
                b_statesFL[5 * A_tmp_tmp + 3] = statesFL[j_tmp + 3];
                b_statesFL[5 * A_tmp_tmp + 4] = 0.0;
            }
            for (A_tmp_tmp = 0; A_tmp_tmp < 5; A_tmp_tmp++) {
                b[static_cast<int>(wf_now + (static_cast<double>(A_tmp_tmp) + 1.0)) -
                  1] = dv1[A_tmp_tmp] - (b_statesFL[A_tmp_tmp] * x1_idx_0 +
                                         b_statesFL[A_tmp_tmp + 5] * x1_idx_1);
            }
        } else {
            j_tmp = b_j << 1;
            for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
                A_tmp = A_tmp_tmp << 2;
                b_A_tmp = 87 * (A_tmp_tmp + j_tmp);
                A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp];
                A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 1];
                A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 2];
                A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 3];
            }
            A[(static_cast<int>(wf_now + 5.0) + 87 * j_tmp) - 1] = 0.0;
            A[(static_cast<int>(wf_now + 5.0) + 87 * (j_tmp + 1)) - 1] = 0.0;
            j_tmp = ((b_j - 1) << 1) + 6;
            for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
                A_tmp = A_tmp_tmp << 2;
                b_A_tmp = 87 * (A_tmp_tmp + j_tmp);
                A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = statesFL[A_tmp];
                A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] = statesFL[A_tmp + 1];
                A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] = statesFL[A_tmp + 2];
                A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] = statesFL[A_tmp + 3];
            }
            A[(static_cast<int>(wf_now + 5.0) + 87 * j_tmp) - 1] = 0.0;
            A[(static_cast<int>(wf_now + 5.0) + 87 * (j_tmp + 1)) - 1] = 0.0;
            for (A_tmp_tmp = 0; A_tmp_tmp < 5; A_tmp_tmp++) {
                A_tmp =
                        static_cast<int>(wf_now + (static_cast<double>(A_tmp_tmp) + 1.0));
                A[(A_tmp + 87 * (b_j + 21)) - 1] = -1.0;
                b[A_tmp - 1] = dv[A_tmp_tmp];
            }
        }
        wf_now += 5.0;
        Q[(b_j + 33 * (b_j + 21)) + 21] = 20000.0;
    }
    //  final state
    for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
        A_tmp = A_tmp_tmp << 2;
        b_A_tmp = 87 * (A_tmp_tmp + 4);
        A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp];
        A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp + 1];
        A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp + 2];
        A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp + 3];
    }
    for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
        A[(static_cast<int>(wf_now + 5.0) + 87 * (A_tmp_tmp + 4)) - 1] = 0.0;
        A_tmp = A_tmp_tmp << 2;
        b_A_tmp = 87 * (A_tmp_tmp + 10);
        A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = statesFL[A_tmp];
        A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] = statesFL[A_tmp + 1];
        A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] = statesFL[A_tmp + 2];
        A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] = statesFL[A_tmp + 3];
    }
    A[static_cast<int>(wf_now + 5.0) + 869] = 0.0;
    A[static_cast<int>(wf_now + 5.0) + 956] = 0.0;
    for (A_tmp_tmp = 0; A_tmp_tmp < 5; A_tmp_tmp++) {
        A_tmp = static_cast<int>(wf_now + (static_cast<double>(A_tmp_tmp) + 1.0));
        A[A_tmp + 2087] = -1.0;
        b[A_tmp - 1] = dv[A_tmp_tmp];
    }
    wf_now += 5.0;
    Q[i + 33 * i] = 20000.0;
    //     %% sigmaFR envelope and slack constrains
    statesFL[0] = 1.6556556940273772;
    statesFL[4] = v;
    statesFL[1] = -1.6556556940273772;
    statesFL[5] = lr;
    statesFL[2] = -1.6556556940273772;
    statesFL[6] = omega_slew;
    statesFL[3] = 1.6556556940273772;
    statesFL[7] = params_idx_1_tmp_tmp;
    i = 25;
    for (b_j = 0; b_j < 3; b_j++) {
        i = b_j + 25;
        if (b_j + 1 == 1) {
            j_tmp = b_j << 1;
            for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
                A_tmp = A_tmp_tmp << 2;
                b_A_tmp = 87 * (A_tmp_tmp + j_tmp);
                A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp];
                A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 1];
                A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 2];
                A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 3];
            }
            A[(static_cast<int>(wf_now + 5.0) + 87 * j_tmp) - 1] = 0.0;
            A[(static_cast<int>(wf_now + 5.0) + 87 * (j_tmp + 1)) - 1] = 0.0;
            for (A_tmp_tmp = 0; A_tmp_tmp < 5; A_tmp_tmp++) {
                A[(static_cast<int>(wf_now + (static_cast<double>(A_tmp_tmp) + 1.0)) +
                   87 * (b_j + 25)) -
                  1] = -1.0;
            }
            for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
                j_tmp = A_tmp_tmp << 2;
                b_statesFL[5 * A_tmp_tmp] = statesFL[j_tmp];
                b_statesFL[5 * A_tmp_tmp + 1] = statesFL[j_tmp + 1];
                b_statesFL[5 * A_tmp_tmp + 2] = statesFL[j_tmp + 2];
                b_statesFL[5 * A_tmp_tmp + 3] = statesFL[j_tmp + 3];
                b_statesFL[5 * A_tmp_tmp + 4] = 0.0;
            }
            for (A_tmp_tmp = 0; A_tmp_tmp < 5; A_tmp_tmp++) {
                b[static_cast<int>(wf_now + (static_cast<double>(A_tmp_tmp) + 1.0)) -
                  1] = dv1[A_tmp_tmp] - (b_statesFL[A_tmp_tmp] * x1_idx_0 +
                                         b_statesFL[A_tmp_tmp + 5] * x1_idx_1);
            }
        } else {
            j_tmp = b_j << 1;
            for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
                A_tmp = A_tmp_tmp << 2;
                b_A_tmp = 87 * (A_tmp_tmp + j_tmp);
                A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp];
                A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 1];
                A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 2];
                A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] =
                        inputsFL_tmp[A_tmp + 3];
            }
            A[(static_cast<int>(wf_now + 5.0) + 87 * j_tmp) - 1] = 0.0;
            A[(static_cast<int>(wf_now + 5.0) + 87 * (j_tmp + 1)) - 1] = 0.0;
            j_tmp = ((b_j - 1) << 1) + 6;
            for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
                A_tmp = A_tmp_tmp << 2;
                b_A_tmp = 87 * (A_tmp_tmp + j_tmp);
                A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = statesFL[A_tmp];
                A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] = statesFL[A_tmp + 1];
                A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] = statesFL[A_tmp + 2];
                A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] = statesFL[A_tmp + 3];
            }
            A[(static_cast<int>(wf_now + 5.0) + 87 * j_tmp) - 1] = 0.0;
            A[(static_cast<int>(wf_now + 5.0) + 87 * (j_tmp + 1)) - 1] = 0.0;
            for (A_tmp_tmp = 0; A_tmp_tmp < 5; A_tmp_tmp++) {
                A_tmp =
                        static_cast<int>(wf_now + (static_cast<double>(A_tmp_tmp) + 1.0));
                A[(A_tmp + 87 * (b_j + 25)) - 1] = -1.0;
                b[A_tmp - 1] = dv[A_tmp_tmp];
            }
        }
        wf_now += 5.0;
        Q[(b_j + 33 * (b_j + 24)) + 24] = 20000.0;
    }
    //  final state
    for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
        A_tmp = A_tmp_tmp << 2;
        b_A_tmp = 87 * (A_tmp_tmp + 4);
        A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp];
        A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp + 1];
        A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp + 2];
        A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] = inputsFL_tmp[A_tmp + 3];
    }
    for (A_tmp_tmp = 0; A_tmp_tmp < 2; A_tmp_tmp++) {
        A[(static_cast<int>(wf_now + 5.0) + 87 * (A_tmp_tmp + 4)) - 1] = 0.0;
        A_tmp = A_tmp_tmp << 2;
        b_A_tmp = 87 * (A_tmp_tmp + 10);
        A[(static_cast<int>(wf_now + 1.0) + b_A_tmp) - 1] = statesFL[A_tmp];
        A[(static_cast<int>(wf_now + 2.0) + b_A_tmp) - 1] = statesFL[A_tmp + 1];
        A[(static_cast<int>(wf_now + 3.0) + b_A_tmp) - 1] = statesFL[A_tmp + 2];
        A[(static_cast<int>(wf_now + 4.0) + b_A_tmp) - 1] = statesFL[A_tmp + 3];
    }
    A[static_cast<int>(wf_now + 5.0) + 869] = 0.0;
    A[static_cast<int>(wf_now + 5.0) + 956] = 0.0;
    for (A_tmp_tmp = 0; A_tmp_tmp < 5; A_tmp_tmp++) {
        A_tmp = static_cast<int>(wf_now + (static_cast<double>(A_tmp_tmp) + 1.0));
        A[A_tmp + 2435] = -1.0;
        b[A_tmp - 1] = dv[A_tmp_tmp];
    }
    wf_now += 5.0;
    Q[i + 33 * i] = 20000.0;
    //     %% tracking objective (for the first steps only)
    //          obj = obj + R*abs(r - u{1}) + R*(r - u{1}).^2;
    //          obj = obj + R*abs(r - u{2}) + R*(r - u{2}).^2;
    for (A_tmp_tmp = 0; A_tmp_tmp < 4; A_tmp_tmp++) {
        i = A_tmp_tmp << 2;
        Q[33 * A_tmp_tmp] = y[i];
        Q[33 * A_tmp_tmp + 1] = y[i + 1];
        Q[33 * A_tmp_tmp + 2] = y[i + 2];
        Q[33 * A_tmp_tmp + 3] = y[i + 3];
    }
    slew[0] = -40.0 * input[5];
    slew[1] = -0.01 * input[6];
    c[0] = slew[0];
    c[29] = 20.0;
    c[1] = slew[1];
    c[30] = 0.005;
    c[2] = slew[0];
    c[31] = 20.0;
    c[3] = slew[1];
    c[32] = 0.005;
    for (i = 0; i < 33; i++) {
        lb[i] = -INFTY;
        ub[i] = INFTY;
    }
    //  abs in obj
    //     %% general dynamics
    //     %% du objective
    for (j = 0; j < 2; j++) {
        params_idx_1_tmp_tmp =
                wf_now + 4.0 * ((static_cast<double>(j) + 1.0) - 1.0);
        b_j = j << 1;
        A[(static_cast<signed char>(params_idx_1_tmp_tmp + 1.0) + 87 * b_j) - 1] =
                -1.0;
        A_tmp = 87 * (b_j + 29);
        A[(static_cast<int>(params_idx_1_tmp_tmp + 1.0) + A_tmp) - 1] = -1.0;
        A[(static_cast<signed char>(params_idx_1_tmp_tmp + 2.0) + 87 * b_j) - 1] =
                1.0;
        A[(static_cast<int>(params_idx_1_tmp_tmp + 2.0) + A_tmp) - 1] = -1.0;
        b[static_cast<signed char>(params_idx_1_tmp_tmp + 1.0) - 1] = -input[5];
        b[static_cast<signed char>(params_idx_1_tmp_tmp + 2.0) - 1] = input[5];
        Aeq[6 * j] = -Bd[b_j];
        i = 6 * (j + 6);
        Aeq[i] = Ad_tmp[b_j];
        params_idx_1_tmp_tmp =
                (wf_now + 4.0 * ((static_cast<double>(j) + 1.0) - 1.0)) + 2.0;
        omega_slew = (wf_now + 4.0 * ((static_cast<double>(j) + 1.0) - 1.0)) + 2.0;
        A_tmp = 87 * (b_j + 1);
        A[(static_cast<signed char>(params_idx_1_tmp_tmp + 1.0) + A_tmp) - 1] =
                -1.0;
        b_A_tmp = 87 * (b_j + 30);
        A[(static_cast<int>(omega_slew + 1.0) + b_A_tmp) - 1] = -1.0;
        A[(static_cast<signed char>(params_idx_1_tmp_tmp + 2.0) + A_tmp) - 1] = 1.0;
        A[(static_cast<int>(omega_slew + 2.0) + b_A_tmp) - 1] = -1.0;
        b[static_cast<signed char>(params_idx_1_tmp_tmp + 1.0) - 1] = -input[6];
        b[static_cast<signed char>(params_idx_1_tmp_tmp + 2.0) - 1] = input[6];
        Aeq[6 * j + 1] = -Bd[b_j + 1];
        Aeq[i + 1] = Ad_tmp[b_j + 1];
        beq[j] = Ad[j] * x1_idx_0 + Ad[j + 2] * x1_idx_1;
        j_tmp = (j + 1) << 1;
        i = j_tmp + 6 * j_tmp;
        Aeq[i] = -Bd[0];
        Aeq[i + 1] = -Bd[1];
        i = j_tmp + 6 * (j_tmp + 1);
        Aeq[i] = -Bd[2];
        Aeq[i + 1] = -Bd[3];
        i = j_tmp + 6 * (b_j + 6);
        Aeq[i] = -(params_idx_0_tmp + 1.0);
        Aeq[i + 1] = -b_params_idx_0_tmp;
        i = j_tmp + 6 * (b_j + 7);
        Aeq[i] = -b_params_idx_1_tmp;
        Aeq[i + 1] = -(params_idx_1_tmp + 1.0);
        i = j_tmp + 6 * (j_tmp + 6);
        Aeq[i] = 1.0;
        Aeq[i + 1] = 0.0;
        i = j_tmp + 6 * (j_tmp + 7);
        Aeq[i] = 0.0;
        Aeq[i + 1] = 1.0;
        c[j] -= b_y[j] * input[j + 7];
    }
    std::memset(&Q_tmp[0], 0, 36U * sizeof(double));
    for (j = 0; j < 6; j++) {
        Q_tmp[j + 6 * j] = dv2[j];
    }
    std::memset(&d[0], 0, 36U * sizeof(double));
    d[12] = -100.0;
    d[19] = -0.3;
    d[26] = -100.0;
    d[33] = -0.3;
    std::memset(&b_d[0], 0, 36U * sizeof(double));
    b_d[2] = -100.0;
    b_d[9] = -0.3;
    b_d[16] = -100.0;
    b_d[23] = -0.3;
    for (A_tmp_tmp = 0; A_tmp_tmp < 36; A_tmp_tmp++) {
        Q_tmp[A_tmp_tmp] = (Q_tmp[A_tmp_tmp] + d[A_tmp_tmp]) + b_d[A_tmp_tmp];
    }
    Q_tmp[35] /= 2.0;
    Q_tmp[28] /= 2.0;
    //     %%
    for (i = 0; i < 6; i++) {
        for (A_tmp_tmp = 0; A_tmp_tmp < 6; A_tmp_tmp++) {
            j_tmp = A_tmp_tmp + 6 * i;
            b_j = A_tmp_tmp + 33 * i;
            params_idx_1_tmp_tmp = Q[b_j] + Q_tmp[j_tmp];
            Q_tmp[j_tmp] = params_idx_1_tmp_tmp;
            Q[b_j] = params_idx_1_tmp_tmp;
        }
        lbA[i] = beq[i];
    }
    for (i = 0; i < 87; i++) {
        lbA[i + 6] = -INFTY;
    }
    for (i = 0; i < 6; i++) {
        ubA[i] = beq[i];
    }
    std::copy(&b[0], &b[87], &ubA[6]);
    for (A_tmp_tmp = 0; A_tmp_tmp < 6; A_tmp_tmp++) {
        for (j_tmp = 0; j_tmp < 33; j_tmp++) {
            Asolver[j_tmp + 33 * A_tmp_tmp] = Aeq[A_tmp_tmp + 6 * j_tmp];
        }
    }
    for (A_tmp_tmp = 0; A_tmp_tmp < 87; A_tmp_tmp++) {
        for (j_tmp = 0; j_tmp < 33; j_tmp++) {
            Asolver[j_tmp + 33 * (A_tmp_tmp + 6)] = A[A_tmp_tmp + 87 * j_tmp];
        }
    }
}


real_t makeStepMPC(real_t input[9],
                 real_t params[9],
                 real_t u[2],
                 SQProblem qrecipeD,
                 real_t r[33],
                 real_t r_dual[33+93],
                 bool firstRun,
                 int infeas_counter[1])
{
    real_t tic = getCPUtime();

    // It's faster to make new matrices than use the previously made one
    real_t Q[1089]; //doesn't change during iterations
    real_t c[33] = {0}; //changes (u_prev)
    real_t Asolver[3069] = {0}; //changes (u_prev + ref + x1)
    real_t lb[33] = {0}; //changes (wr_now)
    real_t ub[33] = {0}; //changes (wr_now)
    real_t lbA[93] = {0}; //changes
    real_t ubA[93] = {0}; //changes

    defineOPC(input, params, Q, c, Asolver, lb, ub, lbA, ubA);

    int_t k = 0;
    int_t pos = 0;
    pos = 0;
    for (k = 0; k < 3069; k++) {
        if (Asolver[k] != 0.0) {
            pos++;
        }
    }
    sparse_int_t A_sparse_rowidx[pos];
    sparse_int_t A_sparse_colidx[33+1];
    real_t A_sparse_d[pos];
    makeSparseFromDense(Asolver, 33, 93, pos, A_sparse_rowidx, A_sparse_colidx, A_sparse_d);

    pos = 0;
    for (k = 0; k < 1089; k++) {
        if (Q[k] != 0.0) {
            pos++;
        }
    }
    sparse_int_t Q_sparse_rowidx[pos];
    sparse_int_t Q_sparse_colidx[33+1];
    real_t Q_sparse_d[pos];
    makeSparseFromDense(Q, 33, 33, pos, Q_sparse_rowidx, Q_sparse_colidx, Q_sparse_d);

    SymDenseMat *Qs_dense = new SymDenseMat(33, 33, 33, Q);
    DenseMatrix *A_dense = new DenseMatrix(93, 33, 33, Asolver);

    SymSparseMat *Qs_sparse = new SymSparseMat(33, 33, Q_sparse_rowidx, Q_sparse_colidx, Q_sparse_d);
    Qs_sparse->createDiagInfo();
    SparseMatrix *A_sparse = new SparseMatrix(93, 33, A_sparse_rowidx, A_sparse_colidx, A_sparse_d);

    real_t my_r[33] = {0}; // for init guess
    my_r[0] = input[7];
    my_r[1] = input[8];
    my_r[2] = input[7];
    my_r[3] = input[8];
    my_r[4] = input[7];
    my_r[5] = input[8];
    my_r[6] = input[0];
    my_r[7] = input[1];
    my_r[8] = input[0];
    my_r[9] = input[1];
    my_r[10] = input[0];
    my_r[11] = input[1];

    qrecipeD.init(Qs_sparse, c, A_sparse, lb, ub, lbA, ubA, nWSR, 0, my_r);

    real_t *out = new real_t[33];
    real_t *out_dual = new real_t[33+93];
    if (!qrecipeD.isSolved())
    {
        u[0] = input[7]; // takes previous step inputs
        u[1] = input[8];
        infeas_counter[0] = infeas_counter[0] + 1;
    }
    else
    {
        qrecipeD.getPrimalSolution(out);
        u[0] = out[0];
        u[1] = out[1];
        for (int i = 0; i < 33; ++i) { // keep for the next iter
            r[i] = out[i];
        }
        qrecipeD.getDualSolution(out_dual);
        for (int i = 0; i < 33+93; ++i) { // keep for the next iter
            r_dual[i] = out_dual[i];
        }
    }

    real_t toc = getCPUtime();

    real_t time = toc-tic;

    return time;
}

void runOneTest(const char* inputFile, const char* paramFile, const char* uFile, const char* tFile)
{
    // data reading for inputs
    std::ifstream data(inputFile);
    std::string line, val;
    std::vector<std::vector<real_t>> input_matrix;
    while (std::getline (data, line)) {
        std::vector<real_t> v;
        std::stringstream s (line);
        while (getline (s, val, ','))
            v.push_back (std::stod (val));
        input_matrix.push_back (v);
    }
    data.close();

    // data reading for params
    std::ifstream data2(paramFile);
    std::vector<std::vector<real_t>> param_matrix;
    while (std::getline (data2, line)) {
        std::vector<real_t> v;
        std::stringstream s(line);
        while (getline(s, val, ','))
            v.push_back(std::stod(val));
        param_matrix.push_back(v);
    }
    data2.close();

    real_t u[2] = {0};
    real_t iter_time = 0;

    std::ofstream u_file;
    u_file.open (uFile);
    std::ofstream t_file;
    t_file.open (tFile);

    SQProblem qrecipeD(33, 93);
    qrecipeD.setPrintLevel(PL_NONE);

    real_t r[33] = {0}; // for init guess
    real_t r_dual[33+93] = {0}; // for init guess

    int infeas_counter[1] = {0};

    for (int_t i = 0; i < param_matrix.size(); ++i) { // param_matrix.size()
        iter_time = makeStepMPC(input_matrix[i].data(), param_matrix[i].data(), u, qrecipeD, r, r_dual, i==0, infeas_counter);
        for (int j = 0; j < 32; ++j) {
            u_file << r[j] << ",";
        }
        u_file << r[32] << std::endl;
        t_file << iter_time << std::endl;
    }
    u_file.close();
    t_file.close();

}

int_t main()
{
    const char * inputBraking = "inputs_beatle_braking.csv";
    const char * paramBraking = "params_beatle_braking.csv";
    const char * uBraking = "u_braking.csv";
    const char * tBraking = "t_braking.csv";
    runOneTest(inputBraking, paramBraking, uBraking, tBraking);
    printf("Braking test is successful \n");

    const char * inputBraking_mu = "inputs_beatle_braking_muSplit.csv";
    const char * paramBraking_mu = "params_beatle_braking_muSplit.csv";
    const char * uBraking_mu = "u_braking_muSplit.csv";
    const char * tBraking_mu = "t_braking_muSplit.csv";
    runOneTest(inputBraking_mu, paramBraking_mu, uBraking_mu, tBraking_mu);
    printf("Braking mu test is successful \n");

    const char * inputAcc = "inputs_acceleration.csv";
    const char * paramAcc = "params_acceleration.csv";
    const char * uAcc = "u_acceleration.csv";
    const char * tAcc = "t_acceleration.csv";
    runOneTest(inputAcc, paramAcc, uAcc, tAcc);
    printf("Acceleration test is successful \n");

    const char * inputSteer = "inputs_steering.csv";
    const char * paramSteer = "params_steering.csv";
    const char * uSteer = "u_steering.csv";
    const char * tSteer = "t_steering.csv";
    runOneTest(inputSteer, paramSteer, uSteer, tSteer);
    printf("Steering test is successful \n");

    return 0;
}
