// SchedulerTest.cpp
#include <scheduler.h>
#include <filehandler.h>
#include <gtest/gtest.h>

TEST(ReadPlan, PlanFileNotExist) {
    vector<Integral> elements;
    bool ret = readPlan(elements, "skdglsdh.txt");

    ASSERT_FALSE(ret);
}

TEST(SchedulerTest, createShells_1) {
    auto origin = vector<double>{2,0,1};
    auto supercell = vector<double>{5,4,3};
    auto unitcell = vector<vector<double>>{
        vector<double>{1,0,0},
        vector<double>{1.1,2,0},
        vector<double>{2.,0,0.5},
    };
    auto r_v = vector<double>{-1,2,3};
    auto r_c = vector<double>{6,-5,4};

    auto shells = createShells(origin,supercell,unitcell,r_v, r_c);

    // number of unitcells
    int N_shouldbe = supercell[0] * supercell[1] * supercell[2];
    int N = 0;
    for (const auto& shell : shells) {
        N = N + shell.size();
        // cout << shell.size() << endl;
    }
    ASSERT_EQ(N, N_shouldbe);

    // max distance to origin has to increase between shells
    double maxDist = 0.0;
    double maxDist_prev = -1.0;
    auto unitcell_T = transpose3x3(unitcell);

    for (const auto& shell : shells) {
        maxDist = 0.0;
        for (const auto& RD : shell) {
            vector<double> vec_shift = matVecMul3x3(unitcell_T, RD);

            double dist = L2Norm(vector<double>{
                r_c[0]-r_v[0]-vec_shift[0],
                r_c[1]-r_v[1]-vec_shift[1],
                r_c[2]-r_v[2]-vec_shift[2]
            });

            // all distances should be larger than the maximal distance of the previous shell
            ASSERT_GT(dist,maxDist_prev);

            maxDist = max(dist,maxDist);
        }
        //cout << maxDist << "\t" << maxDist_prev << endl;
        ASSERT_GT(maxDist, maxDist_prev);
        maxDist_prev = maxDist;
    }

    // check dublicates
    for (size_t i=0; i<shells.size(); i++) {
        for (size_t ii=0; ii<shells[i].size(); ii++) {

            auto RD1 = shells[i][ii];

            for (size_t j=0; j<shells.size(); j++) {
                for (size_t jj=0; jj<shells[j].size(); jj++) {
                    if ((i==j) && (ii==jj)) continue;

                    auto RD2 = shells[j][jj];
                    ASSERT_FALSE( ((RD1[0]==RD2[0]) && (RD1[1]==RD2[1]) && (RD1[2]==RD2[2])) );
                }
            }
        }
    }
}

TEST(SchedulerTest, createShells_2) {
    auto origin = vector<double>{-100,3,4};
    auto supercell = vector<double>{11,12,11};
    auto unitcell = vector<vector<double>>{
        vector<double>{4,0,2.},
        vector<double>{0,2,1},
        vector<double>{-3.,0,5},
    };
    auto r_v = vector<double>{0.3,7,2};
    auto r_c = vector<double>{-15,1,0};

    auto shells{ createShells(origin,supercell,unitcell,r_v, r_c)};

    // number of unitcells
    int N_shouldbe = supercell[0] * supercell[1] * supercell[2];
    //cout << "should be " << N_shouldbe << endl;
    int N = 0;
    for (const auto& shell : shells) {
        N = N + shell.size();
        // cout << shell.size() << endl;
    }

    ASSERT_EQ(N, N_shouldbe);

    // max distance to origin has to increase between shells
    double maxDist = 0.0;
    double maxDist_prev = -1.0;
    auto unitcell_T = transpose3x3(unitcell);

    for (const auto& shell : shells) {
        maxDist = 0.0;
        for (const auto& RD : shell) {
            vector<double> vec_shift = matVecMul3x3(unitcell_T, RD);

            double dist = L2Norm(vector<double>{
                r_c[0]-r_v[0]-vec_shift[0],
                r_c[1]-r_v[1]-vec_shift[1],
                r_c[2]-r_v[2]-vec_shift[2]
            });

            // all distances should be larger than the maximal distance of the previous shell
            ASSERT_GT(dist,maxDist_prev);

            maxDist = max(dist,maxDist);
        }
        // cout << maxDist << "\t" << maxDist_prev << endl;
        ASSERT_GT(maxDist, maxDist_prev);
        maxDist_prev = maxDist;
    }

    // check dublicates
    for (size_t i=0; i<shells.size(); i++) {
        for (size_t ii=0; ii<shells[i].size(); ii++) {

            auto RD1 = shells[i][ii];

            for (size_t j=0; j<shells.size(); j++) {
                for (size_t jj=0; jj<shells[j].size(); jj++) {
                    if ((i==j) && (ii==jj)) continue;

                    auto RD2 = shells[j][jj];

                    ASSERT_FALSE( ((RD1[0]==RD2[0]) && (RD1[1]==RD2[1]) && (RD1[2]==RD2[2])) );

                }
            }
        }
    }
}

vector<Integral> getRandomIntegrals(size_t N) {
    vector<Integral> ret(N);

    for (size_t i=0; i<N; i++) {

        auto indexes = vector<int>(13);
        indexes[0] = i; // to make sure that we dont have dublicates
        for (size_t j=1; j< 4; j++) {
            indexes[j] = rand() % 10; // WF ids must be positive
        }
        for (size_t j=4; j< indexes.size(); j++) {
            indexes[j] = rand() % 20 - 10;
        }
        ret[i] = Integral(indexes);
    }
    return ret;
}


/**
 * @brief Simulates an entire run of the scheduler.
 *
 * It tests whether all integrals are correctly stored.
 * It does not test the correctness of Scheduler::getNextTasks() since it depends on the
 * actual scheduler.
 * Values of the integrals are increasing in steps of 1eV in the order of assignment.
 *
 */
void test_scheduler_basics(Scheduler& myScheduler, size_t num_worker, size_t batch_size, int number_of_prev_known_integrals=0, bool twoTimesSavedInts=false, int N_manual=-1)
{

    ASSERT_FALSE(myScheduler.isFinished());
    ASSERT_EQ(myScheduler.getKnownIntegrals().size(), number_of_prev_known_integrals);

    double value = 0;
    int N = 0;
    map<vector<int>, Integral> knownIntegrals;
    while (! myScheduler.isFinished())
    {
        //cout << "getNextTasks()\n";
        auto tasks = myScheduler.getNextTasks();

        ASSERT_EQ(tasks.size(), num_worker);
        for (size_t worker=0; worker<num_worker; worker++) {
            ASSERT_LE(tasks[worker].size(), batch_size);

            for (size_t i=0; i<tasks[worker].size(); i++) {

                if (tasks[worker][i].isEmpty()) continue;

                tasks[worker][i].value = value;
                tasks[worker][i].setNormal();
                knownIntegrals.insert({tasks[worker][i].indexes, tasks[worker][i]});
                value = value + 1.0;
                N++;
            }
        }

        //cout << "update()\n";
        myScheduler.update(tasks, false);
    }

    auto solved_ints = myScheduler.getKnownIntegrals_asmap();
    auto solved_ints_l = myScheduler.getKnownIntegrals();

    if (N_manual > 0) {
        ASSERT_EQ(solved_ints.size(), N_manual);
    }else if (twoTimesSavedInts) {
        ASSERT_EQ(solved_ints.size(), 2*N + number_of_prev_known_integrals);
    }else{
        ASSERT_EQ(solved_ints.size(), N + number_of_prev_known_integrals);
    }
    ASSERT_EQ(solved_ints.size(), solved_ints_l.size());

    for (const auto& itr : knownIntegrals) {
        auto itr2 = solved_ints.find(itr.first);

        // cout << itr2->second.value << endl;

        ASSERT_NE(itr2, solved_ints.end());
        ASSERT_NEAR(itr.second.value, itr2->second.value, 1e-10);
    }

    for (const auto& i : solved_ints_l) {
        auto itr2 = solved_ints.find(i.indexes);

        ASSERT_NE(itr2, solved_ints.end());
        ASSERT_NEAR(i.value, itr2->second.value, 1e-10);
    }
}

void test_fixedListScheduler(size_t batch_size, size_t num_worker, vector<Integral> plan)
{
    //printPlan(plan, true);
    FixedListScheduler myScheduler(batch_size, num_worker, plan);

    test_scheduler_basics(myScheduler, num_worker, batch_size);

    auto solved_ints_map = myScheduler.getKnownIntegrals_asmap();

    ASSERT_EQ(solved_ints_map.size(), plan.size());

    // check if every element from plan is contained and values are increasing
    double value = 0.0;
    for (const auto& i : plan) {

        //i.print();
        auto itr = solved_ints_map.find(i.indexes);

        ASSERT_TRUE(itr != solved_ints_map.end());
        //cout << itr->second.value << endl;
        //itr->second.print();

        ASSERT_NEAR(itr->second.value, value, 1e-10);
        value = value + 1.0;

    }
}


TEST(SchedulerTest, FixedListScheduler1) {
    auto plan = getRandomIntegrals(rand()%50 + 50);
    test_fixedListScheduler(10, 5, plan);
}


TEST(SchedulerTest, FixedListScheduler2) { // more worker than integrals
    auto plan = getRandomIntegrals(24);
    test_fixedListScheduler(2, 50, plan);
}


TEST(SchedulerTest, FixedListScheduler3) {  // larger batch than integrals
    auto plan = getRandomIntegrals(167);
    test_fixedListScheduler(200, 8, plan);
}

TEST(SchedulerTest, FixedListScheduler4) { // larger batch and workers than integrals
    auto plan = getRandomIntegrals(240);
    test_fixedListScheduler(500, 500, plan);
}

TEST(SchedulerTest, FixedListScheduler5) { // many integrals
    auto plan = getRandomIntegrals(1000000);
    test_fixedListScheduler(10, 12, plan);
}



TEST(SchedulerTest, DensityFixedRScheduler) {

    size_t batch_size = 10;
    size_t num_worker = 4;
    double maxDistance = 35.0;

    auto id_c = vector<int>{3,4,5};
    auto id_v = vector<int>{1,2};


    auto origin = vector<double>{0,3,4};
    auto supercell = vector<double>{11,12,11};
    auto unitcell = vector<vector<double>>{
        vector<double>{4,0,2.},
        vector<double>{0,2,1},
        vector<double>{-3.,0,5},
    };
    auto pos_c = vector<vector<double>>{
        vector<double>{0.3,0.7,2},
        vector<double>{0.3,0.7,2},
        vector<double>{0.3,0.7,2},
    };
    auto pos_v = vector<vector<double>>{
        vector<double>{-1.5,1,0},
        vector<double>{-1.5,1,0},
    };

    auto shells = createShells(origin,supercell,unitcell,pos_v[0], pos_c[0]);

    map<vector<int>, Integral> knownIntegrals;

    Integral a = Integral(1,2,3,4, vector<int>{5,6,7}, vector<int>{8,9,10}, vector<int>{11,12,13});
    a.value = 999.99;
    knownIntegrals.insert({a.indexes, a});

    Integral b = Integral(3,3,1,1, vector<int>{0,3,4}, vector<int>{0,0,0}, vector<int>{0,0,0});
    b.value = 888.88;
    knownIntegrals.insert({b.indexes, b});

    EXPECT_THROW({
        // pos_v and pos_c are exchanged, which should produce an error
         DensityFixedRScheduler(batch_size,num_worker,
            maxDistance, id_c,id_v,shells,pos_v,pos_c,unitcell,knownIntegrals);
    }, runtime_error);


    DensityFixedRScheduler myScheduler(batch_size,num_worker,maxDistance, id_c,id_v,shells,pos_c,pos_v,unitcell,knownIntegrals);

    test_scheduler_basics(myScheduler, num_worker, batch_size,2);

    // check if previous integrals are still contained (not part of the basic test)
    auto solved_ints = myScheduler.getKnownIntegrals_asmap();
    for (const auto& itr : knownIntegrals) {
        auto itr2 = solved_ints.find(itr.first);

        // cout << itr2->second.value << endl;

        ASSERT_NE(itr2, solved_ints.end());
        ASSERT_NEAR(itr.second.value, itr2->second.value, 1e-10);
    }

    // check distance
    auto unitcell_T = transpose3x3(unitcell);
    for (const auto& itr : solved_ints) {
        //itr.second.print();

        if (itr.second.indexes[0] == 1) continue;  // skip previously known integrals, which are outside the sphere

        vector<double> vec_shift = matVecMul3x3(unitcell_T, itr.second.getRD());

        double dist = L2Norm(vector<double>{
            pos_c[0][0]-pos_v[0][0]-vec_shift[0],
            pos_c[0][1]-pos_v[0][1]-vec_shift[1],
            pos_c[0][2]-pos_v[0][2]-vec_shift[2]
        });

        ASSERT_LE(dist, maxDistance);

        // check if new integrals are of type density-density
        ASSERT_EQ(itr.second.indexes[0], itr.second.indexes[1]);
        ASSERT_EQ(itr.second.indexes[2], itr.second.indexes[3]);

        auto Rc = itr.second.getRc();
        ASSERT_EQ(Rc[0], 0);
        ASSERT_EQ(Rc[1], 0);
        ASSERT_EQ(Rc[2], 0);

        auto Rv = itr.second.getRc();
        ASSERT_EQ(Rv[0], 0);
        ASSERT_EQ(Rv[1], 0);
        ASSERT_EQ(Rv[2], 0);

    }

    // TODO check that no integral is forgotten!
}


TEST(SchedulerTest, DensityDensityScheduler) {

    size_t batch_size = 13;
    size_t num_worker = 7;

    auto id_c = vector<int>{3,4,5};
    auto id_v = vector<int>{1,2};


    auto origin = vector<double>{0,3,4};
    auto supercell = vector<double>{11,12,11};
    auto unitcell = vector<vector<double>>{
        vector<double>{4,0,2.},
        vector<double>{0,2,1},
        vector<double>{-3.,0,5},
    };
    auto pos_c = vector<vector<double>>{
        vector<double>{0.3,0.7,2},
        vector<double>{-0.3,0.7,-2},
        vector<double>{0.3,-0.7,2},
    };
    auto pos_v = vector<vector<double>>{
        vector<double>{-1.5,1,0},
        vector<double>{1.5,-1,0},
    };

    auto shells{ createShells(origin,supercell,unitcell,pos_v[0], pos_c[0]) };

    DensityDensityScheduler myScheduler(batch_size,num_worker, id_c,id_v,shells,pos_c,pos_v,unitcell,0.01,999.);

    test_scheduler_basics(myScheduler, num_worker, batch_size);
}

TEST(SchedulerTest, CombindedDensityDensityFixedRScheduler) {

    size_t batch_size = 7;
    size_t num_worker = 13;

    auto id_v = vector<int>{3,4,5};  // id_c and id_v are exchanged compared to DensityDensity test
    auto id_c = vector<int>{1,2};


    auto origin = vector<double>{0,3,4};
    auto supercell = vector<double>{11,12,11};
    auto unitcell = vector<vector<double>>{
        vector<double>{4,0,2.},
        vector<double>{0,2,1},
        vector<double>{-3.,0,5},
    };
    auto pos_v = vector<vector<double>>{
        vector<double>{0.3,0.7,2},
        vector<double>{-0.3,0.7,-2},
        vector<double>{0.3,-0.7,2},
    };
    auto pos_c = vector<vector<double>>{
        vector<double>{-1.5,1,0},
        vector<double>{1.5,-1,0},
    };

    auto shells = createShells(origin,supercell,unitcell,pos_v[0], pos_c[0]);

    CombinedDensityDensityFixedRScheduler myScheduler(batch_size,num_worker, id_c,id_v,shells,pos_c,pos_v,unitcell,0.01, 999.0);

    test_scheduler_basics(myScheduler, num_worker, batch_size);
}


TEST(SchedulerTest, OverlapDensityScheduler) {

    size_t batch_size = 7;
    size_t num_worker = 13;

    double ABSCHARGE_THRESHOLD = 0.1;
    double DISTANCE_THRESHOLD = 100.0;
    double ENERGY_THRESHOLD = 0.01;

    auto origin = vector<double>{0,3,4};
    auto supercell = vector<double>{11,12,11};
    auto unitcell = vector<vector<double>>{
        vector<double>{4,0,2.},
        vector<double>{0,2,1},
        vector<double>{-3.,0,5},
    };

    auto shells{ createShells(origin,supercell,unitcell,vector<double>{0.3,0.7,2}, vector<double>{0.3,0.7,2}) };

    map<Density_descr,Density_indicator> cIndicators{
        { Density_descr(2,3,vector<int>{1,2,3}), Density_indicator(0,1,1, 0.8, -1) },
        { Density_descr(10,6,vector<int>{-1,1,0}), Density_indicator(-2,1,0.4, 0.7, -1) },
        // below ABSCHARGE_THRESHOLD
        { Density_descr(100,60,vector<int>{-1,1,0}), Density_indicator(-2,1,0.4, 0.1-1e-8, -1) }
    };


    map<Density_descr,Density_indicator> vIndicators {
        {Density_descr(3,3,vector<int>{5,1,0-7}), Density_indicator(0,4.5,-0.7, 0.99, -1)}
    };


    OverlapDensityScheduler myScheduler(
        batch_size,num_worker, cIndicators, vIndicators,
        shells, unitcell,
        ABSCHARGE_THRESHOLD, DISTANCE_THRESHOLD, ENERGY_THRESHOLD,
        "OverlapOverlap");

    test_scheduler_basics(myScheduler, num_worker, batch_size, 0, true);

    // TODO: additional tests
}


TEST(SchedulerTest, LocalFieldEffectsScheduler)
{
    size_t batch_size = 7;
    size_t num_worker = 13;
    const double ABSCHARGE_THRESHOLD = 0.1;

    map<Density_descr,Density_indicator> lfe_Indicators {
        {
            Density_descr(2,3,vector<int>{1,2,3}),
            Density_indicator(0.,1.,1., 0.8, -1.)
        },
        {
            Density_descr(10,6,vector<int>{-1,2,0}),
            Density_indicator(-2.,1.,0.4, 0.7, -1.)
        },
        {
            Density_descr(11,7,vector<int>{0,-3,5}),
            Density_indicator(-1.,0.,2.0, 0.86, -1.)
        },
        {
            Density_descr(12,7,vector<int>{0,-3,5}),
            Density_indicator(-1.,0.,2.0, 0.86, -1.)
        },
        {
            Density_descr(11,1,vector<int>{0,0,0}),
            Density_indicator(-1.,0.,2.0, 0.86, -1.)
        },
        {
            Density_descr(1,2,vector<int>{0,-1,0}),
            Density_indicator(-1.,0.,1., 0.6, -1.)
        },
        {  // below ABSCHARGE_THRESHOLD
            Density_descr(100,60,vector<int>{-1,1,0}),
            Density_indicator(-2.,1.,0.4, 0.1-1e-8, -1.)
        }
    };


    LocalFieldEffectsScheduler myScheduler(batch_size, num_worker, lfe_Indicators, ABSCHARGE_THRESHOLD);
    test_scheduler_basics(myScheduler, num_worker, batch_size, 0, true, 36);

    // TODO: additional tests

}