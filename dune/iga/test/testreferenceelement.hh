// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once
#include <execution>
#include <iomanip>
#include <numeric>
#include <random>

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/geometry/referenceelements.hh>
// #include <dune/iga/trimmer/defaulttrimmer/referenceelement.hh>

template <class RE>
auto checkSubEntities(const RE& re) {
  Dune::TestSuite t("checkSubEntities for reference element ");
  auto subEntityName = [](std::size_t codim) -> std::string {
    return codim == 0 ? "element" : codim == 1 ? "edge" : "vertex";
  };
  auto subEntityReport = [&](std::size_t i, std::size_t codim, std::size_t c) {
    std::stringstream ss;

    ss << "The number of " << subEntityName(c) << " of codim " << c << " w.r.t. the " << i << "-th "
       << subEntityName(codim) << " of codim " << c;
    return ss.str();
  };

  for (std::size_t codim = 0; codim <= RE::dimension; ++codim) {
    for (std::size_t i = 0; i < std::size_t(re.size(codim)); ++i) {
      for (std::size_t c = 0; c <= RE::dimension; ++c) {
        auto subEntities = re.subEntities(i, codim, c);
        auto it          = subEntities.begin();
        auto end         = subEntities.end();

        // check size of subEntities range
        auto sizeH = re.size(i, codim, c);

        t.check(sizeH == end - it, "re.size(i, codim, c)== end-it")
            << subEntityReport(i, codim, c) << " is " << sizeH << " but should be " << end - it;
        t.check(std::size_t(sizeH) == subEntities.size(), "std::size_t(re.size(i, codim, c))== subEntities.size()")
            << subEntityReport(i, codim, c) << " is " << sizeH << " but should be " << subEntities.size();

        // check if sub-subentity range is empty if
        // sub-subentity codim c < codim of subentity
        if (c < codim) {
          auto size = re.size(i, codim, c);

          t.check(size == 0, "re.size(i, codim, c)== 0")
              << subEntityReport(i, codim, c) << " is " << size << " but should be " << 0;
        }

        // check if sub-subentity range is singleton if
        // sub-subentity codim c = codim of subentity
        if (c == codim) {
          auto size = re.size(i, codim, c);
          t.check(size == 1, "re.size(i, codim, c)== 1")
              << subEntityReport(i, codim, c) << " is " << size << " but should be " << 1;
        }

        // check if sub-subentity range has multiple
        // entries for sub-subentity codim c > codim of subentity
        if (c > codim) {
          auto size = re.size(i, codim, c);

          t.check(size > 1, "re.size(i, codim, c) > 1")
              << subEntityReport(i, codim, c) << " is " << size << " but should be bigger than " << 1;
        }

        // check if subEntities() is conforming with subEntity()
        for (std::size_t j = 0; j < subEntities.size(); ++j) {
          auto index = std::size_t(re.subEntity(i, codim, j, c));
          t.check(*it == index, "*it== std::size_t(re.subEntity(i, codim, j, c))")
              << "The index of the " << j << "-th " << subEntityName(c) << " of codim " << c << " w.r.t. the " << i
              << "-th " << subEntityName(codim) << " of codim " << codim << " is " << index << " but should be " << *it;
          ++it;
        }

        // check is contains is true for all indices in the range
        for (auto j : subEntities)
          if (not subEntities.contains(j))
            t.check(subEntities.contains(j), "subEntities.contains(j)")
                << "subEntities.contains(" << j << ") returns false but should be true";

        // check if contains is consistent
        std::vector<bool> containedSubEntities(re.size(c), false);
        for (auto j : subEntities)
          containedSubEntities[j] = true;
        for (auto j : Dune::range(std::size_t(re.size(c))))
          t.check(containedSubEntities[j] == subEntities.contains(j),
                  "containedSubEntities[j]==subEntities.contains(j)");

        // after incrementing size times, the range should be exhausted
        t.check(it == end);
      }
    }
  }
  return t;
}

template <class RE>
auto checkCheckInside(const RE& re, int minsamples = 10000, int maxsamples = 100000) {
  Dune::TestSuite t("checkCheckInside for reference element ");

  // int insideCounter=0;

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::default_random_engine random_engine;
  // int samples=1;
  // double oldRatio=-1;
  auto range = Dune::range(maxsamples);
  std::atomic<int> insideCounter{0};
  std::for_each(std::execution::seq, range.begin(), range.end(), [&](int val) {
    double x = unif(random_engine);
    double y = unif(random_engine);
    insideCounter.fetch_add(re.checkInside({x, y}), std::memory_order_relaxed);
  });
  // for (auto i: Dune::range(maxsamples)) {
  //   double x = unif(random_engine );
  //  double y = unif(random_engine );
  //  // std::cout<<"x: "<<x<<" y: "<<y<<std::endl;
  //   if(re.checkInside({x,y}))
  //     ++insideCounter;
  //
  //   if(samples>=minsamples and oldRatio>0 and
  //   std::abs(1-std::abs((static_cast<double>(insideCounter)/samples)/oldRatio))<1e-6 ) {
  //     std::cout<<"Ratio: "<<std::abs(1-std::abs((static_cast<double>(insideCounter)/samples)/oldRatio))<<std::endl;
  //     break;
  //   }
  //   oldRatio=static_cast<double>(insideCounter)/samples;
  //     ++samples;
  // }
  // std::cout<<"insideCounter: "<<insideCounter<<" samples: "<<maxsamples<<std::endl;
  double cubeVolume   = 1.0;
  double insideVolume = cubeVolume * static_cast<double>(insideCounter) / maxsamples;

  t.check(Dune::FloatCmp::eq(re.volume(), insideVolume, 1e-1), "re.volume()==insideVolume")
      << std::setprecision(16) << "re.volume() is " << re.volume() << " but should be " << insideVolume;
  return t;
}

template <typename ReferenceElementImpl, typename ElementTrimDataImpl>
auto checkReferenceElement(const ReferenceElementImpl& refElement, const ElementTrimDataImpl& eleTrimData) {
  Dune::TestSuite t("checkReferenceElement for element ");

  t.subTest(checkSubEntities(refElement));
  t.subTest(checkCheckInside(refElement));
  return t;
}
