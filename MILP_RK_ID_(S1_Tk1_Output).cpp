/*
This code finds all differential characteristics of 13 rounds SKINNY-64 in TK-1 model 
as ∆(S1), ∆(tk11), ∆(output) under the assumption of having a single active bit.
 */ 


#include <fstream>
#include <iostream>
#include <math.h>
#include "gurobi_c++.h" // This library is related to GUROBI solver
using namespace std;

#define ROUNDS (13) // Number of rounds
///////////////////////////////////////////////////////////////////////////////////////////////////
// this is a general MIP callback
class mycallback : public GRBCallback
{
public:
	double lastiter;
	double lastnode;
	int numvars;
	GRBVar* vars;
	ofstream* logfile;
	mycallback(int xnumvars, GRBVar* xvars, ofstream* xlogfile) {
		lastiter = lastnode = -GRB_INFINITY;
		numvars = xnumvars;
		vars = xvars;
		logfile = xlogfile;
	}
protected:
	void callback() {
		try {


			if (where == GRB_CB_MIP) {
				// General MIP callback
				double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
				double objbst = getDoubleInfo(GRB_CB_MIP_OBJBST);
				double objbnd = getDoubleInfo(GRB_CB_MIP_OBJBND);
				int solcnt = getIntInfo(GRB_CB_MIP_SOLCNT);
				if (objbst >= 0) {
					abort();
				}

			}
			else if (where == GRB_CB_MESSAGE) {
				// Message callback
				string msg = getStringInfo(GRB_CB_MSG_STRING);
				*logfile << msg;
			}
		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};
////////////////////////////////////////////////////////////////////////////////////////
/* In this function, the value of all  input/output differentials must be "0". .*/
static void
inOput_to_zero(GRBEnv env, GRBModel& model,
	GRBVar S0[][5][4], GRBVar S1[][5][4], GRBVar S2[][5][4], GRBVar S3[][5][4],
	GRBVar S4[][5][4], GRBVar S5[][5][4], GRBVar S6[][5][4], GRBVar S7[][5][4],
	GRBVar S8[][5][4], GRBVar S9[][5][4], GRBVar S10[][5][4], GRBVar S11[][5][4],
	GRBVar S12[][5][4], GRBVar S13[][5][4], GRBVar S14[][5][4], GRBVar S15[][5][4],
	GRBVar TK0[][5][4], GRBVar TK1[][5][4], GRBVar TK2[][5][4], GRBVar TK3[][5][4],
	GRBVar TK4[][5][4], GRBVar TK5[][5][4], GRBVar TK6[][5][4], GRBVar TK7[][5][4],
	GRBVar TK8[][5][4], GRBVar TK9[][5][4], GRBVar TK10[][5][4], GRBVar TK11[][5][4],
	GRBVar TK12[][5][4], GRBVar TK13[][5][4], GRBVar TK14[][5][4], GRBVar TK15[][5][4],
	GRBConstr m[64], GRBConstr c[64], GRBConstr t[64])
{
	for (int j = 0; j < 4; j++)
	{
		// All differential in S_1 state in 1st round must be zero.
		m[j] = model.addConstr(S0[0][1][j] == 0);
		m[j + 4] = model.addConstr(S1[0][1][j] == 0);
		m[j + 8] = model.addConstr(S2[0][1][j] == 0);
		m[j + 12] = model.addConstr(S3[0][1][j] == 0);
		m[j + 16] = model.addConstr(S4[0][1][j] == 0);
		m[j + 20] = model.addConstr(S5[0][1][j] == 0);
		m[j + 24] = model.addConstr(S6[0][1][j] == 0);
		m[j + 28] = model.addConstr(S7[0][1][j] == 0);
		m[j + 32] = model.addConstr(S8[0][1][j] == 0);
		m[j + 36] = model.addConstr(S9[0][1][j] == 0);
		m[j + 40] = model.addConstr(S10[0][1][j] == 0);
		m[j + 44] = model.addConstr(S11[0][1][j] == 0);
		m[j + 48] = model.addConstr(S12[0][1][j] == 0);
		m[j + 52] = model.addConstr(S13[0][1][j] == 0);
		m[j + 56] = model.addConstr(S14[0][1][j] == 0);
		m[j + 60] = model.addConstr(S15[0][1][j] == 0);
		// All differential in Output state must be zero
		c[j] = model.addConstr(S0[ROUNDS - 1][4][j] == 0);
		c[j + 4] = model.addConstr(S1[ROUNDS - 1][4][j] == 0);
		c[j + 8] = model.addConstr(S2[ROUNDS - 1][4][j] == 0);
		c[j + 12] = model.addConstr(S3[ROUNDS - 1][4][j] == 0);
		c[j + 16] = model.addConstr(S4[ROUNDS - 1][4][j] == 0);
		c[j + 20] = model.addConstr(S5[ROUNDS - 1][4][j] == 0);
		c[j + 24] = model.addConstr(S6[ROUNDS - 1][4][j] == 0);
		c[j + 28] = model.addConstr(S7[ROUNDS - 1][4][j] == 0);
		c[j + 32] = model.addConstr(S8[ROUNDS - 1][4][j] == 0);
		c[j + 36] = model.addConstr(S9[ROUNDS - 1][4][j] == 0);
		c[j + 40] = model.addConstr(S10[ROUNDS - 1][4][j] == 0);
		c[j + 44] = model.addConstr(S11[ROUNDS - 1][4][j] == 0);
		c[j + 48] = model.addConstr(S12[ROUNDS - 1][4][j] == 0);
		c[j + 52] = model.addConstr(S13[ROUNDS - 1][4][j] == 0);
		c[j + 56] = model.addConstr(S14[ROUNDS - 1][4][j] == 0);
		c[j + 60] = model.addConstr(S15[ROUNDS - 1][4][j] == 0);
		// All differential in output of tk_1^1 state must be zero.
		t[j] = model.addConstr(TK0[0][0][j] == 0);
		t[j + 4] = model.addConstr(TK1[0][0][j] == 0);
		t[j + 8] = model.addConstr(TK2[0][0][j] == 0);
		t[j + 12] = model.addConstr(TK3[0][0][j] == 0);
		t[j + 16] = model.addConstr(TK4[0][0][j] == 0);
		t[j + 20] = model.addConstr(TK5[0][0][j] == 0);
		t[j + 24] = model.addConstr(TK6[0][0][j] == 0);
		t[j + 28] = model.addConstr(TK7[0][0][j] == 0);
		t[j + 32] = model.addConstr(TK8[0][0][j] == 0);
		t[j + 36] = model.addConstr(TK9[0][0][j] == 0);
		t[j + 40] = model.addConstr(TK10[0][0][j] == 0);
		t[j + 44] = model.addConstr(TK11[0][0][j] == 0);
		t[j + 48] = model.addConstr(TK12[0][0][j] == 0);
		t[j + 52] = model.addConstr(TK13[0][0][j] == 0);
		t[j + 56] = model.addConstr(TK14[0][0][j] == 0);
		t[j + 60] = model.addConstr(TK15[0][0][j] == 0);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////
/*  This function guarantees the authenticity of all input and output differentials for SKINNY S_box (4*4).
 these constraints are derived from the Algorithm 1 in Sun et.al's paper and SAGE software.*/
static void Sbox(GRBEnv env, GRBModel& model,
	GRBVar x1, GRBVar x2, GRBVar x3, GRBVar x4,
	GRBVar y1, GRBVar y2, GRBVar y3, GRBVar y4)
{
	int mats[24][9] =
	{
		{ -2,3,-3,-2,5,4,1,7,0 },
		{ -1,1,2,1,1,3,0,-2,0 },
		{ 4,3,2,3,-1,-1,-1,-1,0 },
		{ 2,-1,2,0,3,-1,3,-1,0 },
		{ 1,3,1,-2,2,-1,-2,0,2 },
		{ -3,2,1,-2,-1,3,0,1,3 },
		{ -2,-3,2,1,1,-1,3,-1,4 },
		{ 0,1,-2,2,-1,-2,1,1,3 },
		{ 1,-1,-2,-1,1,0,-2,2,4 },
		{ 1,-2,2,-2,-1,1,-1,-2,6 },
		{ 2,1,3,4,-3,2,-1,3,0 },
		{ 0,-1,-1,-1,1,-2,2,-2,5 },
		{ -1,0,-1,-1,-1,1,0,-1,4 },
		{ 1,1,1,2,-2,1,1,0,0 },
		{ 0,-1,-2,-1,1,-1,-2,2,5 },
		{ 0,1,-2,2,-1,-2,-1,-1,5 },
		{ 1,-1,0,-1,1,2,2,2,0 },
		{ 1,0,-1,1,-1,1,0,-1,2 },
		{ -1,-1,-1,1,0,1,0,1,2 },
		{ 3,1,1,-2,2,-1,2,-1,1 },
		{ 0,1,1,0,1,0,-1,0,0 },
		{ -1,-1,1,0,0,-1,1,0,2 },
		{ 1,1,1,0,0,-1,0,0,0 },
		{ 1,0,1,0,-1,-1,-1,0,2 }
	};
	for (int i = 0; i<24; i++)
	{
		model.addConstr((mats[i][0])*x1 + (mats[i][1])*x2 + (mats[i][2])*x3 + (mats[i][3])*x4 +
			(mats[i][4])*y1 + (mats[i][5])*y2 + (mats[i][6])*y3 + (mats[i][7])*y4 >= (-mats[i][8]));
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
/* In this function the input/output variables of SKINNY S-boxes are defined. 
 If we consider S-box state as " Si[j][k][l] " where i=0,...,15 mean a nibble of S-box state,
 "j" is the number of Rounds, "k=0 , 1 " where "k=0" shows  the state before S-box (Input of S-box) and 
 "k=1" shows  the state after S-box (Input of S-box), "l=0,1,2,3" shows the l-bit of "Si" nible.*/
static void sboxx(GRBEnv env, GRBModel& model, int r,
	GRBVar S0[][5][4], GRBVar S1[][5][4], GRBVar S2[][5][4], GRBVar S3[][5][4],
	GRBVar S4[][5][4], GRBVar S5[][5][4], GRBVar S6[][5][4], GRBVar S7[][5][4],
	GRBVar S8[][5][4], GRBVar S9[][5][4], GRBVar S10[][5][4], GRBVar S11[][5][4],
	GRBVar S12[][5][4], GRBVar S13[][5][4], GRBVar S14[][5][4], GRBVar S15[][5][4])
{
	//ROW 1
	Sbox(env, model, S0[r - 1][0][0], S0[r - 1][0][1], S0[r - 1][0][2], S0[r - 1][0][3], S0[r - 1][1][0], S0[r - 1][1][1], S0[r - 1][1][2], S0[r - 1][1][3]);
	Sbox(env, model, S1[r - 1][0][0], S1[r - 1][0][1], S1[r - 1][0][2], S1[r - 1][0][3], S1[r - 1][1][0], S1[r - 1][1][1], S1[r - 1][1][2], S1[r - 1][1][3]);
	Sbox(env, model, S2[r - 1][0][0], S2[r - 1][0][1], S2[r - 1][0][2], S2[r - 1][0][3], S2[r - 1][1][0], S2[r - 1][1][1], S2[r - 1][1][2], S2[r - 1][1][3]);
	Sbox(env, model, S3[r - 1][0][0], S3[r - 1][0][1], S3[r - 1][0][2], S3[r - 1][0][3], S3[r - 1][1][0], S3[r - 1][1][1], S3[r - 1][1][2], S3[r - 1][1][3]);
	//ROW 2
	Sbox(env, model, S4[r - 1][0][0], S4[r - 1][0][1], S4[r - 1][0][2], S4[r - 1][0][3], S4[r - 1][1][0], S4[r - 1][1][1], S4[r - 1][1][2], S4[r - 1][1][3]);
	Sbox(env, model, S5[r - 1][0][0], S5[r - 1][0][1], S5[r - 1][0][2], S5[r - 1][0][3], S5[r - 1][1][0], S5[r - 1][1][1], S5[r - 1][1][2], S5[r - 1][1][3]);
	Sbox(env, model, S6[r - 1][0][0], S6[r - 1][0][1], S6[r - 1][0][2], S6[r - 1][0][3], S6[r - 1][1][0], S6[r - 1][1][1], S6[r - 1][1][2], S6[r - 1][1][3]);
	Sbox(env, model, S7[r - 1][0][0], S7[r - 1][0][1], S7[r - 1][0][2], S7[r - 1][0][3], S7[r - 1][1][0], S7[r - 1][1][1], S7[r - 1][1][2], S7[r - 1][1][3]);
	//ROW 3
	Sbox(env, model, S8[r - 1][0][0], S8[r - 1][0][1], S8[r - 1][0][2], S8[r - 1][0][3], S8[r - 1][1][0], S8[r - 1][1][1], S8[r - 1][1][2], S8[r - 1][1][3]);
	Sbox(env, model, S9[r - 1][0][0], S9[r - 1][0][1], S9[r - 1][0][2], S9[r - 1][0][3], S9[r - 1][1][0], S9[r - 1][1][1], S9[r - 1][1][2], S9[r - 1][1][3]);
	Sbox(env, model, S10[r - 1][0][0], S10[r - 1][0][1], S10[r - 1][0][2], S10[r - 1][0][3], S10[r - 1][1][0], S10[r - 1][1][1], S10[r - 1][1][2], S10[r - 1][1][3]);
	Sbox(env, model, S11[r - 1][0][0], S11[r - 1][0][1], S11[r - 1][0][2], S11[r - 1][0][3], S11[r - 1][1][0], S11[r - 1][1][1], S11[r - 1][1][2], S11[r - 1][1][3]);
	//ROW 4
	Sbox(env, model, S12[r - 1][0][0], S12[r - 1][0][1], S12[r - 1][0][2], S12[r - 1][0][3], S12[r - 1][1][0], S12[r - 1][1][1], S12[r - 1][1][2], S12[r - 1][1][3]);
	Sbox(env, model, S13[r - 1][0][0], S13[r - 1][0][1], S13[r - 1][0][2], S13[r - 1][0][3], S13[r - 1][1][0], S13[r - 1][1][1], S13[r - 1][1][2], S13[r - 1][1][3]);
	Sbox(env, model, S14[r - 1][0][0], S14[r - 1][0][1], S14[r - 1][0][2], S14[r - 1][0][3], S14[r - 1][1][0], S14[r - 1][1][1], S14[r - 1][1][2], S14[r - 1][1][3]);
	Sbox(env, model, S15[r - 1][0][0], S15[r - 1][0][1], S15[r - 1][0][2], S15[r - 1][0][3], S15[r - 1][1][0], S15[r - 1][1][1], S15[r - 1][1][2], S15[r - 1][1][3]);

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// In this function the structure of tweakey in TK-1 model is modeled. 
static void TK01(GRBEnv env, GRBModel& model, int r,
	GRBVar TK0[][5][4], GRBVar TK1[][5][4], GRBVar TK2[][5][4], GRBVar TK3[][5][4],
	GRBVar TK4[][5][4], GRBVar TK5[][5][4], GRBVar TK6[][5][4], GRBVar TK7[][5][4],
	GRBVar TK8[][5][4], GRBVar TK9[][5][4], GRBVar TK10[][5][4], GRBVar TK11[][5][4],
	GRBVar TK12[][5][4], GRBVar TK13[][5][4], GRBVar TK14[][5][4], GRBVar TK15[][5][4])
{
	for (int j = 0; j<4; j++)
	{
		model.addConstr(TK0[r][0][j] == TK8[r + 1][0][j]);
		model.addConstr(TK1[r][0][j] == TK9[r + 1][0][j]);
		model.addConstr(TK2[r][0][j] == TK10[r + 1][0][j]);
		model.addConstr(TK3[r][0][j] == TK11[r + 1][0][j]);
		model.addConstr(TK4[r][0][j] == TK12[r + 1][0][j]);
		model.addConstr(TK5[r][0][j] == TK13[r + 1][0][j]);
		model.addConstr(TK6[r][0][j] == TK14[r + 1][0][j]);
		model.addConstr(TK7[r][0][j] == TK15[r + 1][0][j]);
		model.addConstr(TK8[r][0][j] == TK2[r + 1][0][j]);
		model.addConstr(TK9[r][0][j] == TK0[r + 1][0][j]);
		model.addConstr(TK10[r][0][j] == TK4[r + 1][0][j]);
		model.addConstr(TK11[r][0][j] == TK7[r + 1][0][j]);
		model.addConstr(TK12[r][0][j] == TK6[r + 1][0][j]);
		model.addConstr(TK13[r][0][j] == TK3[r + 1][0][j]);
		model.addConstr(TK14[r][0][j] == TK5[r + 1][0][j]);
		model.addConstr(TK15[r][0][j] == TK1[r + 1][0][j]);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// In this function the XOR operetion is modeled. 
static void XORR(GRBEnv env, GRBModel& model, int num_rounds,
	GRBVar A0[][5][4], GRBVar A1[][5][4], GRBVar A2[][5][4], GRBVar A3[][5][4], int a,
	GRBVar B0[][5][4], GRBVar B1[][5][4], GRBVar B2[][5][4], GRBVar B3[][5][4], int b,
	GRBVar C0[][5][4], GRBVar C1[][5][4], GRBVar C2[][5][4], GRBVar C3[][5][4], int c)
{
	GRBVar D[1][4][4];

	for (int i = 0; i<4; i++)
		for (int j = 0; j<4; j++)
			D[0][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	model.update();


	for (int j = 0; j<4; j++)
	{
		model.addConstr(A0[num_rounds][a][j] + B0[num_rounds][b][j] + C0[num_rounds][c][j] == 2 * D[0][0][j]);
		model.addConstr(A1[num_rounds][a][j] + B1[num_rounds][b][j] + C1[num_rounds][c][j] == 2 * D[0][1][j]);
		model.addConstr(A2[num_rounds][a][j] + B2[num_rounds][b][j] + C2[num_rounds][c][j] == 2 * D[0][2][j]);
		model.addConstr(A3[num_rounds][a][j] + B3[num_rounds][b][j] + C3[num_rounds][c][j] == 2 * D[0][3][j]);
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	// Open log file
	ofstream logfile("cb.log");
	if (!logfile.is_open()) {
		cout << "Cannot open cb.log for callback message" << endl;
		return 1;
	}
	GRBEnv *env = 0;
	GRBVar *vars = 0;

	//using fo print the input/output of characteristics
	ofstream outputFile("output.dat", ios::out);
	if (!outputFile)
	{
		cerr << "some thing wrong during opening file!" << endl;
		exit(1);
	}

	try
	{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		// Set variables
		GRBVar S0[ROUNDS][5][4];
		GRBVar S1[ROUNDS][5][4];
		GRBVar S2[ROUNDS][5][4];
		GRBVar S3[ROUNDS][5][4];
		GRBVar S4[ROUNDS][5][4];
		GRBVar S5[ROUNDS][5][4];
		GRBVar S6[ROUNDS][5][4];
		GRBVar S7[ROUNDS][5][4];
		GRBVar S8[ROUNDS][5][4];
		GRBVar S9[ROUNDS][5][4];
		GRBVar S10[ROUNDS][5][4];
		GRBVar S11[ROUNDS][5][4];
		GRBVar S12[ROUNDS][5][4];
		GRBVar S13[ROUNDS][5][4];
		GRBVar S14[ROUNDS][5][4];
		GRBVar S15[ROUNDS][5][4];
		GRBVar TK0[ROUNDS][5][4];
		GRBVar TK1[ROUNDS][5][4];
		GRBVar TK2[ROUNDS][5][4];
		GRBVar TK3[ROUNDS][5][4];
		GRBVar TK4[ROUNDS][5][4];
		GRBVar TK5[ROUNDS][5][4];
		GRBVar TK6[ROUNDS][5][4];
		GRBVar TK7[ROUNDS][5][4];
		GRBVar TK8[ROUNDS][5][4];
		GRBVar TK9[ROUNDS][5][4];
		GRBVar TK10[ROUNDS][5][4];
		GRBVar TK11[ROUNDS][5][4];
		GRBVar TK12[ROUNDS][5][4];
		GRBVar TK13[ROUNDS][5][4];
		GRBVar TK14[ROUNDS][5][4];
		GRBVar TK15[ROUNDS][5][4];
		GRBConstr m[64];
		GRBConstr c[64];
		GRBConstr t[64];

		for (int i = 0; i < ROUNDS; i++) {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 4; k++)
				{
					S0[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S1[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S2[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S3[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S4[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S5[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S6[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S7[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S8[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S9[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S10[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S11[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S12[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S13[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S14[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					S15[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK0[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK1[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK2[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK3[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK4[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK5[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK6[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK7[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK8[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK9[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK10[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK11[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK12[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK13[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK14[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					TK15[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				}
			}
		}
		model.update();
		//Add S-Boxes to MILP model.
		for (int i = 1; i <= ROUNDS; i++)
			sboxx(env, model, i, S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15);
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Add tweakeys to the rounds in TK1 model
		for (int i = 0; i < ROUNDS; i++)
		{
			XORR(env, model, i, S0, S1, S2, S3, 1, TK0, TK1, TK2, TK3, 0, S0, S1, S2, S3, 2);
			XORR(env, model, i, S4, S5, S6, S7, 1, TK4, TK5, TK6, TK7, 0, S4, S5, S6, S7, 2);
		}
		for (int i = 0; i < ROUNDS; i++)
			for (int j = 0; j < 4; j++)
			{
				model.addConstr(S8[i][1][j] == S8[i][2][j]);
				model.addConstr(S9[i][1][j] == S9[i][2][j]);
				model.addConstr(S10[i][1][j] == S10[i][2][j]);
				model.addConstr(S11[i][1][j] == S11[i][2][j]);
				model.addConstr(S12[i][1][j] == S12[i][2][j]);
				model.addConstr(S13[i][1][j] == S13[i][2][j]);
				model.addConstr(S14[i][1][j] == S14[i][2][j]);
				model.addConstr(S15[i][1][j] == S15[i][2][j]);
			}
		// Add P operetion in the tweakey to the MILP model.
		if (ROUNDS >= 2)
		{
			for (int i = 0; i < ROUNDS - 1; i++)
				TK01(env, model, i, TK0, TK1, TK2, TK3, TK4, TK5, TK6, TK7, TK8, TK9, TK10, TK11, TK12, TK13, TK14, TK15);
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////
//shiftrow operation
		for (int i = 0; i < ROUNDS; i++)
			for (int j = 0; j < 4; j++)
			{
				model.addConstr(S0[i][2][j] == S0[i][3][j]);
				model.addConstr(S1[i][2][j] == S1[i][3][j]);
				model.addConstr(S2[i][2][j] == S2[i][3][j]);
				model.addConstr(S3[i][2][j] == S3[i][3][j]);
				model.addConstr(S4[i][2][j] == S5[i][3][j]);
				model.addConstr(S5[i][2][j] == S6[i][3][j]);
				model.addConstr(S6[i][2][j] == S7[i][3][j]);
				model.addConstr(S7[i][2][j] == S4[i][3][j]);
				model.addConstr(S8[i][2][j] == S10[i][3][j]);
				model.addConstr(S9[i][2][j] == S11[i][3][j]);
				model.addConstr(S10[i][2][j] == S8[i][3][j]);
				model.addConstr(S11[i][2][j] == S9[i][3][j]);
				model.addConstr(S12[i][2][j] == S15[i][3][j]);
				model.addConstr(S13[i][2][j] == S12[i][3][j]);
				model.addConstr(S14[i][2][j] == S13[i][3][j]);
				model.addConstr(S15[i][2][j] == S14[i][3][j]);
//////////////////////////////////////////////////////////////////////////////////////////
//mixcolumns  operation
				model.addConstr(S0[i][3][j] == S4[i][4][j]);
				model.addConstr(S1[i][3][j] == S5[i][4][j]);
				model.addConstr(S2[i][3][j] == S6[i][4][j]);
				model.addConstr(S3[i][3][j] == S7[i][4][j]);
			}
		for (int i = 0; i < ROUNDS; i++)
		{
			XORR(env, model, i, S4, S5, S6, S7, 3, S8, S9, S10, S11, 3, S8, S9, S10, S11, 4);
			XORR(env, model, i, S0, S1, S2, S3, 3, S8, S9, S10, S11, 3, S12, S13, S14, S15, 4);
			XORR(env, model, i, S12, S13, S14, S15, 4, S12, S13, S14, S15, 3, S0, S1, S2, S3, 4);
		}
// After second round, the output of the r-th round must equal to the input of (r+1)-th round
		if (ROUNDS >= 2)
			for (int i = 0; i < ROUNDS - 1; i++)
				for (int j = 0; j < 4; j++)
				{
					model.addConstr(S0[i][4][j] == S0[i + 1][0][j]);
					model.addConstr(S1[i][4][j] == S1[i + 1][0][j]);
					model.addConstr(S2[i][4][j] == S2[i + 1][0][j]);
					model.addConstr(S3[i][4][j] == S3[i + 1][0][j]);
					model.addConstr(S4[i][4][j] == S4[i + 1][0][j]);
					model.addConstr(S5[i][4][j] == S5[i + 1][0][j]);
					model.addConstr(S6[i][4][j] == S6[i + 1][0][j]);
					model.addConstr(S7[i][4][j] == S7[i + 1][0][j]);
					model.addConstr(S8[i][4][j] == S8[i + 1][0][j]);
					model.addConstr(S9[i][4][j] == S9[i + 1][0][j]);
					model.addConstr(S10[i][4][j] == S10[i + 1][0][j]);
					model.addConstr(S11[i][4][j] == S11[i + 1][0][j]);
					model.addConstr(S12[i][4][j] == S12[i + 1][0][j]);
					model.addConstr(S13[i][4][j] == S13[i + 1][0][j]);
					model.addConstr(S14[i][4][j] == S14[i + 1][0][j]);
					model.addConstr(S15[i][4][j] == S15[i + 1][0][j]);
				}
////////////////////////////////////////////////////////////////////////////////////////////////////////
//At first, we consider all input/output differentials zero 
		inOput_to_zero(env, model, S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, TK0, TK1, TK2, TK3, TK4, TK5, TK6, TK7, TK8, TK9, TK10, TK11, TK12, TK13, TK14, TK15, m, c, t);
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Define the Objective Function
		GRBLinExpr OBJS, OBJTK1, OBJ1S, OBJ1TK1, OBJ1TK2, OBJ1TK3, OBJnS, OBJnTK1, OBJu;
		model.update();
		for (int i = 0; i < ROUNDS; i++)
			for (int j = 0; j < 4; j++)
				OBJS += (S0[i][0][j] + S1[i][0][j] + S2[i][0][j] + S3[i][0][j] + S4[i][0][j] + S5[i][0][j] + S6[i][0][j] + S7[i][0][j] + S8[i][0][j] + S9[i][0][j] + S10[i][0][j] + S11[i][0][j] + S12[i][0][j] + S13[i][0][j] + S14[i][0][j] + S15[i][0][j]);
		for (int j = 0; j < 4; j++) {
			OBJ1S += (S0[0][0][j] + S1[0][0][j] + S2[0][0][j] + S3[0][0][j] + S4[0][0][j] + S5[0][0][j] + S6[0][0][j] + S7[0][0][j] + S8[0][0][j] + S9[0][0][j] + S10[0][0][j] + S11[0][0][j] + S12[0][0][j] + S13[0][0][j] + S14[0][0][j] + S15[0][0][j]);
			OBJ1TK1 += (TK0[0][0][j] + TK1[0][0][j] + TK2[0][0][j] + TK3[0][0][j] + TK4[0][0][j] + TK5[0][0][j] + TK6[0][0][j] + TK7[0][0][j] + TK8[0][0][j] + TK9[0][0][j] + TK10[0][0][j] + TK11[0][0][j] + TK12[0][0][j] + TK13[0][0][j] + TK14[0][0][j] + TK15[0][0][j]);
		}
// At least one variable must be nonzero
		model.addConstr(OBJ1S + OBJ1TK1 >= 1);

		model.update();
		model.setObjective(OBJS, GRB_MINIMIZE);
////////////////////////////////////////////////////////////////////////////////////////////////////////
// model.optimize();
/* In the following, we optimize MILP model with considering all inputs/outputs as zero
and then set them one by one the value of "1" and optimize the model. for each input/utput
that the model is infeasible, we have an impossible differential characteristic with that
input/utput differential.*/
		int numvars = model.get(GRB_IntAttr_NumVars);
		vars = model.getVars();
		mycallback cb = mycallback(numvars, vars, &logfile);
		model.setCallback(&cb);
		model.optimize();
		std::cout << "------------------------------------------------------------------------------------\n ";
		for (int i = 0; i < 63; i++)
		{
			//We set the "i"-th bit of S_1^1 state to "1".
			m[i].set(GRB_DoubleAttr_RHS, 1.0);
			for (int j = 0; j < 63; j++)
			{
				//We set the "j"-th bit of tk_1^1 state to "1".
				t[j].set(GRB_DoubleAttr_RHS, 1.0);
				for (int k = 0; k < 63; k++)
				{
					//We set the "k"-th bit of Output state to "1".
					c[k].set(GRB_DoubleAttr_RHS, 1.0);


					model.update();
					// Solve the model
					model.optimize();
					if (model.get(GRB_IntAttr_Status) == 3)
					{
						std::cout << "model is infeasible" << "\n";
						//Print the result if the (i,j,k) is Im_RK characteristic.
						outputFile << " The difference of (" << i<<","<<j<<","<<k << ")-th bit of (Input, tk_1^1, Output) is IM_RK for "
							<<ROUNDS <<"-round of SKINNY in TK1 model.\n";
					}
					//We set the "k"-th bit of Output state to "0".
					c[k].set(GRB_DoubleAttr_RHS, 0.0);
					model.update();
				}
				//We set the "j"-th bit of tk_1^1 state to "0".
				t[j].set(GRB_DoubleAttr_RHS, 0.0);
				model.update();
			}
			//We set the "k"-th bit of Output state to "0".
			m[i].set(GRB_DoubleAttr_RHS, 0.0);
			model.update();
		}
	}
//////////////////////////////////////////////////////////////////////////////////////////
	catch (GRBException e)
	{
		std::cout << "Error code = "
			<< e.getErrorCode()
			<< std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...)
	{
		std::cout << "Exception during optimization"
			<< std::endl;
	}

	return 0;
}