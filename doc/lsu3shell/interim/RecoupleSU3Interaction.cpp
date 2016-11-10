#include<SU3ME/SU3InteractionRecoupler.h>

using namespace std;

int main(int argc, char* argv[])
{
	
	if (argc != 3)
	{
		cout << "Usage: RecoupleSU3NNInteraction <input filename> <output filename prefix>" << endl;

		cout << "This program takes as an input list of SU(3) tensors operators in form" << endl;
		cout << "n1 n2 n3 n4 IRf Sf IRi Si IR0 S0" << endl;
		cout << "{a_{pp} 	a_{nn} 		a_{pn}}_{rho0=0}" << endl;
		cout <<	"{a_{pp}	a_{nn}		a_{pn}}_{rho0=1}" << endl;
		cout << "				. " << endl;
		cout << " 				. " << endl;
		cout << "				. " << endl;
		cout << "{a_{pp}	a_{nn}		a_{pn}}_{rho0=rho0max}" << endl;
		cout << "and represents the following tensor:" << endl;
		cout << "[ {a+_n1 x a+_n2}^{IRf}Sf x {a_n3 x a_n4}^{IRi}Si ]^{IR0} S0" << endl << endl; 

		cout << "Name of the file with SU3 tensors is required as an input argument." << endl;
		return 0;
	}

	recoupler_nonscalar::SU3InteractionRecoupler Recoupler;
	string sInputFileName(argv[1]);
	string sOutputFileName(argv[2]);
	std::ifstream inter_file(sInputFileName.c_str());
	size_t nTensors2Recouple;

	if (!inter_file)
	{
		cerr << "Could not open the file '" << sInputFileName << "'" << endl;
		return 0;
	}
	cout << "reading tensors from the input file '" << sInputFileName << "'" << endl;

	InitSqrtLogFactTables();

	while (true)
	{
		if (inter_file.eof())
		{
			break;
		}

		SU3::LABELS IRf, IRi, Irrep0;
		size_t nCoeffs, index;
		int Sf, Si, S0;
		int n1, n2, n3, n4;
		int k0, k0max, L0, rho0max, rho0;

		inter_file >> n1;
		inter_file >> n2;
		inter_file >> n3;
		inter_file >> n4;

		inter_file >> IRf;
		inter_file >> Sf;
		inter_file >> IRi;
		inter_file >> Si;
		inter_file >> Irrep0;
		inter_file >> S0;

		if (inter_file.eof())
		{
			break;
		}

		rho0max = Irrep0.rho;
		L0 = S0;

//		k0max = SU3::kmax(Irrep0, L0/2);   // HERE WAS AN ERROR: one needs to provide L0/2 instead of L0 !!!
		k0max = 1;
		nCoeffs = 3*k0max*rho0max;
		std::vector<double> CoeffsPPNNPN(nCoeffs, 0.0);
		for (size_t i = 0; i < nCoeffs; ++i)
		{
			inter_file >> CoeffsPPNNPN[i];
		}

		if (std::count_if(CoeffsPPNNPN.begin(), CoeffsPPNNPN.end(), Negligible) == nCoeffs)
		{
			cerr << "All coefficients of the tensor [" << n1 << " " << n2 << " " << n3 << " " << n4 << " " << IRf << " Sf=" << Sf << "/2 " << IRi << " Si=" << Si << "/2 " << Irrep0;
			cerr << " S0=" << S0 << "/2 ] are equal to zero!" << endl;
			exit(EXIT_FAILURE);
		}
///////////////////////////////////////////////////////////////////////////////////////////
		cout << "Recoupling: " << endl;
		cout << "n1 = " << n1 << " n2 = " << n2 << " n3 = " << n3 << " n4 = " << n4;
		cout << "\t" << IRf << " x " << IRi << " -> " << Irrep0;
		cout << "\t" << "Sf = " << Sf << "/2  Si = " << Si << "/2 S0 = " << S0 << "/2" << endl;
		for (k0 = 0, index = 0; k0 < k0max; ++k0)
		{	
			for (rho0 = 0; rho0 < rho0max; ++rho0, index += 3)
			{
				cout << CoeffsPPNNPN[index + PP] << " " << CoeffsPPNNPN[index + NN] << " " << CoeffsPPNNPN[index + PN] << endl;
			}
		}
///////////////////////////////////////////////////////////////////////////////////////////

		char n1n2n3n4[4];
		n1n2n3n4[0] = n1; n1n2n3n4[1] = n2; n1n2n3n4[2] = n3; n1n2n3n4[3] = n4;

//	Note that InsertTTY1 uses formulas for recoupling which do not depend on J0
//	and M0 (see A.35-A.43) and hence these quantum labeks do not need to be specified. 
		int bSuccess = Recoupler.Insert_adad_aa_Tensor(n1n2n3n4, SU3xSU2::LABELS(IRf, Sf), SU3xSU2::LABELS(IRi, Si), SU3xSU2::LABELS(Irrep0, S0), CoeffsPPNNPN);
		if (!bSuccess)
		{
			cout << "IMPLEMENT: ";
		
			cout << "n1 = " << n1 << " n2 = " << n2 << " n3 = " << n3 << " n4 = " << n4;
			cout << "\t" << IRf << " x " << IRi << " -> " << Irrep0;
			cout << "\t" << "Sf = " << Sf << "/2  Si = " << Si << "/2 S0 = " << S0 << "/2";
			cout << "\tk0 = " << k0 << " L0 = " << L0 << "/2" << endl;
			for (k0 = 0, index = 0; k0 < k0max; ++k0)
			{	
				for (rho0 = 0; rho0 < rho0max; ++rho0, index += 3)
				{
					cout << CoeffsPPNNPN[index + PP] << " " << CoeffsPPNNPN[index + NN] << " " << CoeffsPPNNPN[index + PN] << endl;
				}
			}
			cout << endl;
			return 0;
		}
	}

	Recoupler.RemoveTensorsWithAmplitudesLess(9.0e-6); // remove tensors with amplitude lees than 1.0e-5
	Recoupler.PrintInfo();

	string sOutputFileNamePPNN(sOutputFileName);
	string sOutputFileNamePN(sOutputFileName);
	sOutputFileNamePN += ".PN";
	sOutputFileNamePPNN += ".PPNN";
	Recoupler.Save(sOutputFileNamePPNN, TENSOR::PPNN);
	Recoupler.Save(sOutputFileNamePN, TENSOR::PN);

	Recoupler.Show2();
}
