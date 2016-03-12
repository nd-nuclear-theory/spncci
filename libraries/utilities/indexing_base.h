/****************************************************************
  indexing_base.h

  Defines template base classes for indexing quantum mechanical states
  arranged into subspaces (defined by good symmetry quantum numbers).

  General scheme:

    The foundational class is the "subspace" type, representing a
    subspace with good quantum numbers (e.g., LSJT).  This class sets
    up and stores indexing of state quantum numbers within that
    subspace.

    Access to this indexing is via instantiating a "state".  A state
    has to refer to a "subspace" in order to know its quantum numbers.
    Its definition therefore depends upon the "subspace" type.  The
    actual information content of a state consists of a pointer to the
    subspace and an index representing the state within that subspace.

    Going in the other direction in terms of the hierarchy of
    mathematical structures, a "space" is a collection "subspaces".
    Its definition therefore *also* depends upon the "subspace" type.

    Finally, in order to keep track of the matrix elements of an
    operator (of good JP quantum numbers), it is necessary to
    enumerate the allowed sectors, i.e., the pairs of "subspaces"
    within the "space" which are connected under the angular momentum
    and parity selection rules.  These sectors are indexed by a
    "sector enumeration".  Its definition depends upon the "space"
    type, and therefore indirectly on the "subspace" type.

    So notice that the hiererchy of dependencies of template
    definitions 

        state <- subspace -> space -> sectors

    does not quite parallel the mathematical hiererchy

        state -> subspace -> space -> sectors.

  Possible extensions:

    Arguably, each subspace and space should have defined a
    TruncationLabelsType, with a corresponding data member and
    accessors.  The constructor would then record the truncation
    scheme for posterity and later reference.

    For standardized debugging output, add a String conversion to
    interface of each type (state, subspace, ...), then use this in
    standardized Print functions.


  Language: C++11
                                 
  Mark A. Caprio, University of Notre Dame.

  11/21/15 (mac): Created, abstracted from code in indexing_lsjt.
  3/5/16 (mac): Update header comment.
  3/9/16 (mac): Relax protection of BaseSpace indexing containers 
    from private to protected to support sp3rlib.
  3/11/16 (mac): Add subspace size() method and space 
    TotalDimension() method.

****************************************************************/

#ifndef INDEXING_BASE_H_
#define INDEXING_BASE_H_

#include <cassert>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

namespace shell {


  ////////////////////////////////////////////////////////////////
  // generic subspace
  ////////////////////////////////////////////////////////////////

  // BaseSubspace -- holds indexing of states sharing common subspace
  // quantum numbers
  //
  // The derived class is expected to set up a constructor and
  // friendlier accessors for the individual labels.
  //
  // Template arguments:
  //   tSubspaceLabelsType : tuple for subspace labels, e.g., std::tuple<int,int,int,int,int>
  //   tStateLabelsType : tuple for state labels, e.g., std::tuple<int> 
  
  template <typename tSubspaceLabelsType, typename tStateLabelsType> 
    class BaseSubspace
  {

  public:

    ////////////////////////////////////////////////////////////////
    //  common type definitions
    ////////////////////////////////////////////////////////////////

    typedef tSubspaceLabelsType SubspaceLabelsType; 
    typedef tStateLabelsType StateLabelsType; 

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    // default constructor
    //   Implicitly invoked by derived class.
  BaseSubspace() : dimension_(0) {}

    // copy constructor -- synthesized


    ////////////////////////////////////////////////////////////////
    // retrieval
    ////////////////////////////////////////////////////////////////

    const SubspaceLabelsType& GetSubspaceLabels() const
    // const SubspaceLabels& GetSubspaceLabels() const
    {
      return labels_;
    }

    const StateLabelsType& GetStateLabels(int index) const
    // const StateLabels& GetStateLabels(int index) const
    {
      return state_table_[index];
    }

    int LookUpStateIndex(const StateLabelsType& state_labels) const
      {
	return lookup_.at(state_labels);
      };

    int size() const
    // DEPRECATED in favor of size()
    {
      return dimension_;
    }

    int Dimension() const
    // DEPRECATED in favor of size()
    {
      return dimension_;
    }

  protected:

    ////////////////////////////////////////////////////////////////
    // state label push (for initial construction)
    ////////////////////////////////////////////////////////////////

    void PushStateLabels(const StateLabelsType& state_labels)
    {
      lookup_[state_labels] = dimension_; // index for lookup
      state_table_.push_back(state_labels);  // save state
      dimension_++;
    };
    
    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // subspace properties
    SubspaceLabelsType labels_;
    int dimension_;

    // state labels (accessible by index)
    std::vector<StateLabelsType> state_table_;

    // state index lookup by labels
    std::map<StateLabelsType,int> lookup_;

  };

  ////////////////////////////////////////////////////////////////
  // generic state realized within subspace
  ////////////////////////////////////////////////////////////////

  // BaseState -- realization of a state withinin a given subspace
  //
  // The derived class is expected to set up a constructor and
  // friendlier accessors for the individual labels.
  //
  //
  // The space (and the indexing it provides) is *not* copied into the
  // state but rather stored by pointer reference.  It should
  // therefore exist for the lifetime of the state object.
  //
  // Template arguments:
  //   tSubspaceType : subspace type in which this state lives
  
  template <typename tSubspaceType> 
    class BaseState
  {

  public:

    ////////////////////////////////////////////////////////////////
    // common typedefs constructors
    ////////////////////////////////////////////////////////////////

    typedef tSubspaceType SubspaceType; 

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    // default constructor -- disabled
    BaseState();

    // copy constructor -- synthesized

    // constructors

    BaseState(const SubspaceType& subspace, int index)
      // Construct state by index.
      : subspace_ptr_(&subspace), index_(index)
      {
	assert(ValidIndex());
      }

    BaseState(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
      {

	// debugging: Delegation to BaseState(subspace,index)
	// fails. Argument index is saved in index_ as far as the
	// subordinate BaseState(...,index) call is concerned, but
	// after return to present calling constructor, index_
	// contains garbage.  (Why???)

	// BaseState(subspace,index);

	subspace_ptr_ = &subspace;
	index_ = subspace.LookUpStateIndex(state_labels);

      }

    ////////////////////////////////////////////////////////////////
    // retrieval
    ////////////////////////////////////////////////////////////////

    const SubspaceType& Subspace() const {return *subspace_ptr_;}  // subspace in which state lies
    const typename SubspaceType::StateLabelsType& GetStateLabels() const 
    {
      return Subspace().GetStateLabels(index());
    }

    int index() const {return index_;}

    ////////////////////////////////////////////////////////////////
    // generic iteration support -- disabled
    //
    // currently DISABLED, pending decision about whether or not
    // this is a good thing
    //
    // Example:
    //
    //   for (RelativeStateLSJT state(space); state.ValidIndex(); ++state)
    //     {
    //   	std::cout << state.index() << " " << state.N() << std::endl;
    //     };
    //
    ////////////////////////////////////////////////////////////////

    BaseState(const SubspaceType& subspace);
    // Constructs state, defaulting to 0th state in space.
    // Meant for use with iterator-style iteration over states.
    //  {
    //	space_ptr = subspace;
    //	index_ = 0;
    //  }

    BaseState& operator ++();
    // Provides prefix increment operator.
    //
    // Meant for use with iterator-style iteration over states.
    // {
    // 	++index_;
    // 	return *this;
    // }


  private:

    ////////////////////////////////////////////////////////////////
    // validation
    ////////////////////////////////////////////////////////////////

    int ValidIndex() const
    // Verifies whether or not state indexing lies within allowed
    // dimension.
    //
    // For use on construction or possibly with iteration.
    {
      return index() < Subspace().size();
    }

  private:

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // const SubspaceType& subspace_;  // subspace in which state lies
    const SubspaceType* subspace_ptr_;  // subspace in which state lies
    int index_;   // 0-based index within space

  };



  ////////////////////////////////////////////////////////////////
  // generic space
  ////////////////////////////////////////////////////////////////

  // BaseSpace -- container to hold subspaces of type S with reverse
  // lookup by subspace labels
  //
  // Template arguments:
  //   tSubspaceType (typename) : type for subspace

  template <typename tSubspaceType>
    class BaseSpace
    {
    public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      typedef tSubspaceType SubspaceType;
      
      ////////////////////////////////////////////////////////////////
      // subspace lookup and retrieval
      ////////////////////////////////////////////////////////////////

      int LookUpSubspaceIndex(const typename SubspaceType::SubspaceLabelsType& subspace_labels) const
      {
	
	// diagnostic output for failed lookup
	if (!lookup_.count(subspace_labels))
	  {

	    // std::cout << "Space lookup: label not found " << subspace_labels << std::endl;
	    std::cout << "Map size: " << lookup_.size() << std::endl;
	    // for (auto& elem : lookup_) // doesn't work
	    //for (auto iter =lookup_.begin(); iter != lookup_.end(); ++iter)
	    //  {
	    //	const auto& ikey = (*iter).first;
	    //	std::cout << std::get<0>(ikey) << "," << std::get<1>(ikey) << "," << std::get<2>(ikey) << "," << std::get<3>(ikey);
	    //	std::cout << " : ";
	    //	std::cout << (*iter).second;
	    //	std::cout << " match ";
	    //	std::cout << (subspace_labels == ikey);
	    //	std::cout << std::endl;
	    //  }
	    std::cout << std::endl;
	  }

	return lookup_.at(subspace_labels);
      };

      const SubspaceType& LookUpSubspace(const typename SubspaceType::SubspaceLabelsType& labels) const
      {
	return subspaces_[LookUpSubspaceIndex(labels)];
      };

      const SubspaceType& GetSubspace(int i) const
      {
	return subspaces_[i];
      };

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

      int size() const
      {
	return subspaces_.size();
      };

      int TotalDimension() const
      // Sum dimensions of subspaces.
      {
        int dimension = 0;
        for (int i=0; i<size(); ++i)
          dimension += GetSubspace(i).size();
      }

    protected:

      ////////////////////////////////////////////////////////////////
      // subspace push (for initial construction)
      ////////////////////////////////////////////////////////////////

      void PushSubspace(const SubspaceType& subspace)
      {
	lookup_[subspace.GetSubspaceLabels()] = subspaces_.size(); // index for lookup
	subspaces_.push_back(subspace);  // save space
      };

      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      // subspaces (accessible by index)
      std::vector<SubspaceType> subspaces_;

      // subspace index lookup by labels
      std::map<typename SubspaceType::SubspaceLabelsType,int> lookup_;
      
    };

  ////////////////////////////////////////////////////////////////
  // generic sector storage
  ////////////////////////////////////////////////////////////////

  // Sector -- sector index pair
  //
  // Derived from pair<int,int> so that can be easily used as sortable key.
  // Public inheritance so comparisons work.

  class Sector
    : public std::pair<int,int> 
  {

  public:

    // constructor
  Sector(int index2, int index1) 
    : pair(index2,index1) {}

    // accessors (aliases for pair members)

    int index2() const {return first;}
    int index1() const {return second;}

  };

  // BaseSectors -- container to hold sectors (as pairs of
  // indices) with reverse lookup by space labels
  //
  // Effectively an invertible function ZxZ->Z.
  //
  // It is not clear there is any actual need for this class to be
  // templatized, since tSpaceType is not actually used, but it
  // maintains logical parallelism to the other base classes in this
  // header.
  //
  // Template arguments:
  //   tSpaceType (typename) : type for space

  template <typename tSpaceType>
    class BaseSectors
  {
  public:

    ////////////////////////////////////////////////////////////////
    // common typedefs
    ////////////////////////////////////////////////////////////////

    typedef tSpaceType SpaceType;

    // constructor to be provided by derived class

    ////////////////////////////////////////////////////////////////
    // subspace lookup and retrieval
    ////////////////////////////////////////////////////////////////

    const Sector& GetSector(int i) const
    {
      return sectors_[i];
    };

    int LookUpSectorIndex(const Sector& sector) const
    {
      return lookup_.at(sector);
    };

    ////////////////////////////////////////////////////////////////
    // size retrieval
    ////////////////////////////////////////////////////////////////

    int size() const
    {
      return sectors_.size();
    };

  protected:

    ////////////////////////////////////////////////////////////////
    // sector push (for initial construction)
    ////////////////////////////////////////////////////////////////
    

    void PushSector(const Sector& sector)
    {
      lookup_[sector] = sectors_.size(); // index for lookup
      sectors_.push_back(sector);  // save sector
    };

    ////////////////////////////////////////////////////////////////
    // internal storage
    ////////////////////////////////////////////////////////////////
    
    // sectors (accessible by index)
    std::vector<Sector> sectors_;

    // sector index lookup by subspace indices
    std::map<Sector,int> lookup_;
  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
