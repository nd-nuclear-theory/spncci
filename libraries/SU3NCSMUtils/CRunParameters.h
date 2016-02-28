#ifndef RUN_PARAMETERS_H
#define RUN_PARAMETERS_H
#include <SU3ME/InteractionPPNN.h>
#include <SU3ME/proton_neutron_ncsmSU3Basis.h>
#include <SU3ME/CInteractionPN.h>

class CRunParameters
{
	private:
	friend class boost::serialization::access;
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & hw_;
        ar & ncsmModelSpace_;
		// I do not serialize data with input interactions/operators
		// since only root process needs to know that.
    }
	private:
	proton_neutron::ModelSpace ncsmModelSpace_;
	double hw_;
	std::vector<std::pair<std::string, float> > two_body_PPNN_;
	std::vector<std::pair<std::string, float> > two_body_PN_;
	std::vector<std::pair<std::string, float> > one_body_;
	private:
	inline double A() const {return (double)(Z() + N());}

	inline double TrelCoeff() const {return (2.0*hw_)/A();}
	inline double NcmCoeff() const {return 1.0/A();}
	double VcoulCoeff() const { const static double fmass = 938.92; return sqrt(fmass/(double)938.093)*sqrt(hw_/(double)(10.0)); }

	inline void AddOneBodyOperator(const std::string& observable_1b, const double dcoeff = 1.0) {one_body_.push_back(std::make_pair(observable_1b, dcoeff));}
	inline void AddTwoBodyOperatorPPNN(const std::string& observable_2b_PPNN, const double dcoeff = 1.0) {two_body_PPNN_.push_back(std::make_pair(observable_2b_PPNN, dcoeff));}
	inline void AddTwoBodyOperatorPN(const std::string& observable_2b_PN, const double dcoeff = 1.0) {two_body_PN_.push_back(std::make_pair(observable_2b_PN, dcoeff));}

	void AddTrel();
	void AddVcoul();
	void AddVnn();
	void AddNcm(const double lambda);
	void AddAB00(const double lambda);

	public:
	void LoadRunParameters(const char* run_params_file_name);
	public:
	inline double hw() const {return hw_;}
	inline int N() const {return ncsmModelSpace_.number_of_neutrons();}
	inline int Z() const {return ncsmModelSpace_.number_of_protons();}
	inline int Nmax() const {return ncsmModelSpace_.back().N();}
	inline const proton_neutron::ModelSpace& GetModelSpace() const {return ncsmModelSpace_;}
	public:
	void LoadInteractionTerms(int my_rank, CInteractionPPNN& interactionPPNN, CInteractionPN& interactionPN);
//	new interface:	
	void LoadHamiltonian(int my_rank, const std::string& hamiltonian_file_name, CInteractionPPNN& interactionPPNN, CInteractionPN& interactionPN, int A);

	inline double TrelCoeff(int hw, int A) const {return (2.0*hw)/A;}
	void AddTrel(int hw, int A);

	void AddVcoul(int hw);
	inline double VcoulCoeff(int hw) const {double fmass = 938.92; return sqrt(fmass/(double)938.093)*sqrt(hw/(double)(10.0));}
	
	void AddNcm(const double lambda, int A);
	inline double NcmCoeff(int A) const {return 1.0/A;}
};

#endif
