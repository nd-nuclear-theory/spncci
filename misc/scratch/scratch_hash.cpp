    inline friend std::size_t hash_value(const OperatorLabelsU3ST& v)
    {
      // std::size_t seed = 0;
      // boost::hash_combine(seed,v.N0_);
      // boost::hash_combine(seed,v.x0_);
      // boost::hash_combine(seed,v.S0_);
      // boost::hash_combine(seed,v.T0_);
      // boost::hash_combine(seed,v.g0_);
      // return seed;

      boost::hash<OperatorLabelsU3ST::KeyType> hasher;
      return hasher(v.Key());

    }
