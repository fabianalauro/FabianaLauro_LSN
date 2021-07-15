//
//
//
//
//  Created by Fabiana on 06/07/21.
//

#include "Objects.hpp"


//City ---------------------------------------------------------------------------
City::City()
{
    m_label = 0 ;
    m_first = "false" ;
    m_x = 0. ;
    m_y = 0. ;
}


City::City (int label, bool first, double x, double y)
{
    m_label = label ;
    m_first = first ;
    m_x = x ;
    m_y = y ;
}

City& City::operator= (const City& C)
{
    m_label = C.label() ;
    m_first = C.first() ;
    m_x = C.x() ;
    m_y = C.y() ;
    return *this ;
}

void City::SetX (double x)
{
    m_x = x;
}
    
void City::SetY (double y)
{
    m_y = y ;
}

void City::SetL (int l) 
{
    m_label = l ;
}


void City::PrintCity () 
{
    std::cout << "CITY " << m_label << ": fist=" << m_first << " x=" << m_x << " y=" << m_y << std::endl ;
}

void City::PrintCity (std::ofstream &out) 
{
    out << m_label << "\t" << m_x << "\t" << m_y << std::endl ;
}



//Trip ---------------------------------------------------------------------------

Trip::Trip()
{
    m_ncities = 0. ;
    m_length = 0.;
    m_rnd = NULL ;
    m_cities = std::vector<City> (m_ncities) ;
}


Trip::Trip (int N, Random * gen )
{
    m_ncities = N ;
    m_length = 0.;
    m_rnd = gen ;
    m_cities = std::vector<City> (m_ncities) ;
}


void Trip::SetCities (std::vector<City> V)
{
    m_cities = V ;
    return;
}

void Trip::SetPositions (double * x, double * y) 
{
    for (int i = 0; i < m_ncities; i++)
    {
        m_cities [i].SetX(x[i]);
        m_cities [i].SetY(y[i]);
    }
}

void Trip::SetPositions (int * l, double * x, double * y) 
{
    for (int i = 0; i < m_ncities; i++)
    {
        m_cities [i].SetL(l[i]);
        m_cities [i].SetX(x[i]);
        m_cities [i].SetY(y[i]);
    }
}


Trip& Trip::operator= (const Trip& T)
{
    m_ncities = T.ncities() ;
    m_length = T.length() ;
    m_rnd = T.rnd() ;
    m_cities.clear() ;
    for (int i=0; i<T.ncities(); i++)
        m_cities.push_back(T.cities()[i]);
    
    return *this ;
}

bool Trip::CheckFirst () 
{
    return m_cities[0].first() ;
}

void Trip::AddCity (City C, int index) 
{
    if (index >= m_ncities)
        std::cerr << "WARNING: the index of the City you're trying to add is greater than ncities" << std::endl ;
    m_cities[index] = C ;
}

void Trip::measure ()
{
    double dist = 0. ;
    for (int i = 0; i < m_ncities-1; i++)
    {
        dist += sqrt ( pow (m_cities[i].x() - m_cities[i+1].x(), 2) + pow (m_cities[i].y() - m_cities[i+1].y(), 2) )   ;
    }
    dist += sqrt ( pow (m_cities[m_ncities-1].x() - m_cities[0].x(), 2) + pow (m_cities[m_ncities-1].y() - m_cities[0].y(), 2) ) ;
    
    m_length = dist ;
}

int Trip::Pbc (int index, int size)
{
    if (index < size && index > 0)
        return index;
    else if (index < 0)
    {
        while (index < 0)
            index += (size-1) ;
        index++;
        return index;

    }
    else
    {
        while (index >= size)
            index -= (size-1) ;
        return index;

    }
}

void Trip::Print () 
{
    for (int i=0; i<m_cities.size(); i++)
        m_cities[i].PrintCity() ;
}

void Trip::Print (std::ofstream &out) 
{
    for (int i=0; i<m_cities.size(); i++)
    {
        out << i << "\t" ;
        m_cities[i].PrintCity(out) ;
    }
}

void Trip::Swap () 
{
    
    int I1 = trunc (m_rnd->Rannyu(1, m_ncities));
    int I2 = trunc (m_rnd->Rannyu(1, m_ncities));
    City foo = m_cities [Pbc(I1, m_ncities)] ;
    m_cities [Pbc(I1, m_ncities)] = m_cities [Pbc(I2, m_ncities)];
    m_cities[Pbc(I2, m_ncities)] = foo ;
}

void Trip::SwapCont() 
{
    int index = trunc (m_rnd->Rannyu(1, m_ncities));
    City foo = m_cities [index] ;
    m_cities [index] = m_cities [Pbc(index+1, m_ncities)] ;
    m_cities [Pbc(index+1, m_ncities)] = foo ;
}

void Trip::SwapM() 
{
    int m = trunc ( m_rnd->Rannyu(2, m_ncities/2.) ); //Number of elements to swap
    int I1 = trunc ( m_rnd->Rannyu(1, m_ncities-m) ) ; 
    int I2 = trunc ( m_rnd-> Rannyu(m+I1, m_ncities+(I1-1)-m) ) ;
    std::vector<City> foo ;
    for (int i=0; i<m; i++)
    {
        foo.push_back (m_cities[I1+i]) ;
        m_cities[I1+i] = m_cities[Pbc(I2+i, m_ncities)] ;
        m_cities[Pbc(I2+i, m_ncities)] = foo[i] ;
    }
}

void Trip::Shift() 
{
    int I = trunc ( m_rnd->Rannyu(1, m_ncities) ); //Index of the first to shift
    int m = trunc ( m_rnd->Rannyu(1, m_ncities-I) ); //Number of elements to shift
    int n = trunc ( m_rnd->Rannyu(1, m_ncities-1) ); //Number of positions to shift

    std::vector<City>::iterator It ;
    std::vector<City> appo ;
    for (int i=0; i<m; i++)
        appo.push_back(m_cities[Pbc(I+i, m_ncities)]);  

    It = m_cities.begin() ;
    m_cities.erase(It+I, It+I+m);

    It = m_cities.begin() ;
    m_cities.insert(It+Pbc(I+n, m_cities.size()), appo.begin(), appo.end()) ;

}

void Trip::Reverse()
{
    int I = trunc ( m_rnd->Rannyu(1, m_ncities) ); //Index of the first reverse
    int m = trunc ( m_rnd->Rannyu(2, m_ncities) ); //Number of elements to reverse

    std::vector<City>::iterator It = m_cities.begin() + I ;
    std::vector<City>::reverse_iterator rIt ;
    std::vector<City> appo ;
    for (int i=0; i<m; i++)
        appo.push_back(m_cities[Pbc(I+i, m_ncities)]); 
    for (rIt = appo.rbegin(); rIt < appo.rend(); rIt++)
    {
        *It = *rIt ;
        It ++ ;
        if (It == m_cities.end())
            It = m_cities.begin() + 1;
    }
}

//Annealing ---------------------------------------------------------------------------

Annealing::Annealing () 
{
    Trip m_TheTrip ; 
    Trip m_NewTrip ; 
    m_rnd = NULL ; 
    m_beta = 0. ; 
    m_count = 0. ; 
    m_alpha = 0. ; 
    Trip m_best ; 
}


Annealing::Annealing (int ncities, Trip T, Random * gen) 
{
    m_TheTrip = T ; 
    Trip m_NewTrip (ncities, gen); 
    m_rnd = gen ; 
    m_beta = 0. ; 
    m_count = 0. ; 
    m_alpha = 0. ; 
    m_best = T; 
}


void Annealing::Mutate () 
{
    m_NewTrip = m_TheTrip ;
    double r = m_rnd->Rannyu(0,5) ;
    if (trunc(r) < 1)
        m_NewTrip.Swap () ;
    else if (r > 1 && r < 2)
        m_NewTrip.SwapCont () ;
    else if (r > 2 && r < 3)
        m_NewTrip.SwapM () ;
    else if (r > 3 && r < 4)
        m_NewTrip.Shift () ;
    else if (r > 4 && r < 5)
        m_NewTrip.Reverse () ;
}


void Annealing::MRT2 () 
{   
    m_count ++ ;
    this->Mutate ();
    m_NewTrip.measure();
    m_TheTrip.measure();
    double Boltzmann = exp ( -m_beta * (m_NewTrip.length() - m_TheTrip.length()) );
    m_alpha = fmin (1, Boltzmann) ;
    if (m_rnd->Rannyu() < m_alpha)
        m_TheTrip = m_NewTrip ;
    if (m_TheTrip.length() < m_best.length() )
    {
        m_best = m_TheTrip ;
        m_bestcount = m_count ;
    }
}

void Annealing::Print (bool The) 
{
    if (The == true)
        m_TheTrip.Print();
    else 
        m_NewTrip.Print();
}

void Annealing::Print (std::ofstream &out) 
{
    m_TheTrip.Print (out) ;
    out << 32 << "\t";
    m_TheTrip.cities()[0].PrintCity(out);
}

void Annealing::PrintEvoBest (std::ofstream &out) 
{
    out << m_count << "\t" << m_beta << "\t" << m_best.length() << std::endl  ;
}

void Annealing::PrintEvo (std::ofstream &out) 
{
    out << m_count << "\t" << m_beta << "\t" << m_TheTrip.length() << "\t" << m_bestcount << std::endl  ;
}


//Population ---------------------------------------------------------------------------

Population::Population ()
{
    m_ntrips = 0;
    m_rnd = NULL ;
    m_generation = 0 ;
    m_bestgen = 0 ;
    m_trips = std::vector<Trip> (m_ntrips);
    m_deliveryroom = std::vector<Trip> (m_ntrips);
    m_migr = std::vector<int> (0) ;
}

Population::Population (int N, Random * gen)
{
    m_ntrips = N;
    m_rnd = gen ;
    m_generation = 0 ;
    m_bestgen = 0 ;
    m_trips = std::vector<Trip> (m_ntrips);
    m_deliveryroom = std::vector<Trip> (m_ntrips);
    m_migr = std::vector<int> (0) ;

}

Population::Population (int N, Random * gen, int ncores) 
{
    m_ntrips = N;
    m_rnd = gen ;
    m_generation = 0 ;
    m_bestgen = 0 ;
    m_trips = std::vector<Trip> (m_ntrips);
    m_deliveryroom = std::vector<Trip> (m_ntrips);
    m_migr = std::vector<int> (ncores) ;
    for (int i=0; i<ncores; i++)
        m_migr[i] = i ;
    
}

void Population::SwapMigr () 
{
    int I1 = trunc (m_rnd->Rannyu(0, m_migr.size()));
    int I2 = trunc (m_rnd->Rannyu(0, m_migr.size()));
    int foo = m_migr[I1] ;
    m_migr[I1] = m_migr[I2];
    m_migr[I2] = foo ;
}

void Population::SetMigr (int * v) 
{
    for (int i=0; i<m_migr.size() ; i++)
        m_migr[i] = v[i];
}

void Population::Swap (int index1, int index2, bool pop) 
{
    if (index1 == index2)
        return;
    else
    {
        if (pop == true)
        {
            Trip foo (m_trips[index1].ncities(), m_rnd) ;
            foo = m_trips[index1];
            m_trips[index1] = m_trips[index2] ;
            m_trips[index2] = foo ;
        }
        else
        {
            Trip foo (m_deliveryroom[index1].ncities(), m_rnd) ;
            foo = m_deliveryroom[index1];
            m_deliveryroom[index1] = m_deliveryroom[index2] ;
            m_deliveryroom[index2] = foo ;
        }
   }
}

void Population::AddTrip (Trip T, int index) 
{
    if (index >= m_ntrips)
        std::cerr << "WARNING: the index of the Trip you're trying to add is greater than ntrips" << std::endl ;
    m_trips[index] = T ;
}

void Population::Sort (bool pop) 
{
    int start = 0 ;
    int lowest ;
    while (start < m_ntrips-1)
    {
        lowest = start ;
        for (int i=start; i<m_ntrips; i++)
        {
            if (pop == true)
            {
                m_trips[i].measure() ;
                if ( m_trips[i].length() < m_trips[lowest].length() )
                {
                    lowest = i ;
                }
            }
            else
            {
                m_deliveryroom[i].measure() ;
                if ( m_deliveryroom[i].length() < m_deliveryroom[lowest].length() )
                {
                    lowest = i ;
                } 
            }
        }
        Swap (start, lowest, pop) ;
        start ++ ;
    }
    if (pop == true)
    {
        if (m_trips[0].length() < m_best.length())
        {
            m_best = m_trips[0];
            m_bestgen = m_generation ;
        }
    }

}


void Population::PrintLength (bool pop) 
{
    for (int i=0; i<m_ntrips; i++)
    {
        if (pop == true)
        {
            m_trips[i].measure();
            std::cout << "TRIP #" << i << "   length=" << m_trips[i].length() << std::endl;
        }
        else
        {
            m_deliveryroom[i].measure();
            std::cout << "TRIP #" << i << "   length=" << m_deliveryroom[i].length() << std::endl;
        }
    }
}

void Population::PrintBest()
{
    std::cout << "GENERATION #" << m_generation << "\t BEST TRIP: d=" << m_best.length() << "\t BEST GEN: " << m_bestgen << std::endl;
}

void Population::PrintBest(std::ofstream &out)
{
    out << m_generation << "\t" << m_best.length() << "\t\t" << m_bestgen << std::endl;
}

void Population::PrintSol (std::ofstream &out)
{
    m_trips[0].Print(out) ;
    out << 32 << "\t";
    m_trips[0].cities()[0].PrintCity(out);
}


void Population::PrintSol ()
{
    m_trips[0].Print() ;
    std::cout << 32 << "\t";
    m_trips[0].cities()[0].PrintCity();

}




void Population::PrintAve (std::ofstream &out) 
{
    this->Sort(true) ;
    double ave = 0. ;
    for (int i=0; i<m_ntrips/2; i++)
    {
        ave += m_trips[i].length();
    }
    ave /= m_ntrips/2. ;
    out << m_generation << "\t" << ave << std::endl;

}

void Population::Print () 
{
    for (int i=0; i<m_ntrips; i++)
    {
        std::cout << "TRIP #" << i << std::endl;
        m_trips[i].Print() ;
    }
}

void Population::Print (int n) 
{
    for (int i=0; i<n; i++)
    {
        std::cout << "TRIP #" << i << std::endl;
        m_trips[i].Print() ;
    }
}

void Population::Selection (double p)
{
    for (int i=0; i<m_ntrips; i+=2)
    {
        int I1 = trunc ( m_ntrips * pow(m_rnd->Rannyu(), p) ) ;
        int I2 = 0 ;
        do
            I2 = trunc ( m_ntrips * pow(m_rnd->Rannyu(), p) ) ;
        while (I2 == I1);

        m_deliveryroom[i] = m_trips[I1] ;
        m_deliveryroom[i+1] = m_trips[I2];
    }
    
}

void Population::Crossover (Trip &P1, Trip &P2)
{
    bool found ;
    int i = P1.ncities()/2;
    int j = 0 ;
    while (i<P1.ncities() )
    {   
        found = false;
        for (int k=0; k<P1.ncities()/2; k++) 
        {
            if (P2.cities()[j].label() == P1.cities()[k].label())
            {
                found = true ;
                break ;
            }
        }
        if (found == false)
        {
            P1.AddCity(P2.cities()[j], i) ;
            i++;
        }
        j++ ;
        if (j > P2.ncities())
            std::cerr << "attention! j>ncities!!!!!!!! I=" << i << std::endl;
    }

}

void Population::GiveBirth (double p) 
{
    if (m_rnd->Rannyu() < p)
    {
        Trip foo ( m_deliveryroom[0].ncities(), m_deliveryroom[0].rnd() ) ;
        for (int i=0 ; i<m_ntrips; i+=2)
        {
            foo = m_deliveryroom[i] ;
            this->Crossover( m_deliveryroom[i], m_deliveryroom[i+1] );
            this->Crossover( m_deliveryroom[i+1], foo );
        }
    }
}

void Population::Mutation(double p)
{
    for (int i=0; i<m_ntrips; i++)
    {
        if (m_rnd->Rannyu() < p)
            m_deliveryroom[i].SwapCont ();
        if (m_rnd->Rannyu() < p)
            m_deliveryroom[i].Swap () ;
        if (m_rnd->Rannyu() < p)
            m_deliveryroom[i].SwapM () ;
        if (m_rnd->Rannyu() < p)
            m_deliveryroom[i].Shift ();
        if (m_rnd->Rannyu() < p)
            m_deliveryroom[i].Reverse ();
    }
}

void Population::NewGen ()
{
    this->Sort(true);
    this->Sort(false) ;
    for (int i=0; i<m_ntrips/4; i++)
    {
        m_trips[m_ntrips/4+i] = m_deliveryroom[i] ;
    }
    m_generation ++;
}

void Population::Migration (int rank)
{
    int ncities = m_trips[0].ncities();
    int s_l [ncities];
    int r_l [ncities];
    double s_x [ncities];
    double s_y [ncities];
    double r_x [ncities];
    double r_y [ncities];

    for (int i=0; i<ncities; i++)
    {
        s_l[i] = m_trips[0].cities()[i].label();
        s_x[i] = m_trips[0].cities()[i].x();
        s_y[i] = m_trips[0].cities()[i].y();
    }
    /*
    if (rank == m_migr[1])
    {
        std::cout << "rank " << m_migr[1] << " sending" << std::endl ;
        std::cout << "x = " ;
        for (int i=0; i<ncities; i++)
            std:: cout << s_x[i] << " " ;
        std::cout << std::endl ;
        std::cout << "y = " ;
        for (int i=0; i<ncities; i++)
            std:: cout << s_y[i] << " " ;
        std::cout << std::endl ;
    }
*/
    int send13x = 1;
    int recv13x = 2;
    int send02x = 3;
    int recv02x = 4;
    int send13y = 5;
    int recv13y = 6;
    int send02y = 7;
    int recv02y = 8;
    int send13l = 9 ;
    int send02l = 10 ;
    int recv13l = 11 ;
    int recv02l = 12 ;
    
    MPI_Status s13x, s31x, s02x, s20x, s13y, s31y, s02y, s20y, s13l, s31l, s20l, s02l;

    if(rank == m_migr[0])
    {
        std::cout <<"From " << rank << " to "  << m_migr[2] << std::endl;
        MPI_Send(&s_x[0], ncities, MPI_DOUBLE, m_migr[2], send02x, MPI_COMM_WORLD);
        MPI_Recv(&r_x[0], ncities, MPI_DOUBLE, m_migr[2], recv02x, MPI_COMM_WORLD, &s20x);
        MPI_Send(&s_y[0], ncities, MPI_DOUBLE, m_migr[2], send02y, MPI_COMM_WORLD);
        MPI_Recv(&r_y[0], ncities, MPI_DOUBLE, m_migr[2], recv02y, MPI_COMM_WORLD, &s20y);
        MPI_Send(&s_l[0], ncities, MPI_INTEGER, m_migr[2], send02l, MPI_COMM_WORLD);
        MPI_Recv(&r_l[0], ncities, MPI_INTEGER, m_migr[2], recv02l, MPI_COMM_WORLD, &s20l);
    }

    else if (rank == m_migr[1])
    {
        std::cout <<"From " << rank << " to "  << m_migr[3] << std::endl;
        MPI_Send(&s_x[0], ncities, MPI_DOUBLE, m_migr[3], send13x, MPI_COMM_WORLD);
        MPI_Recv(&r_x[0], ncities, MPI_DOUBLE, m_migr[3], recv13x, MPI_COMM_WORLD, &s13x);
        MPI_Send(&s_y[0], ncities, MPI_DOUBLE, m_migr[3], send13y, MPI_COMM_WORLD);
        MPI_Recv(&r_y[0], ncities, MPI_DOUBLE, m_migr[3], recv13y, MPI_COMM_WORLD, &s13y);
        MPI_Send(&s_l[0], ncities, MPI_INTEGER, m_migr[3], send13l, MPI_COMM_WORLD);
        MPI_Recv(&r_l[0], ncities, MPI_INTEGER, m_migr[3], recv13l, MPI_COMM_WORLD, &s13l);

    }

    else if (rank == m_migr[2])
    {
        std::cout <<"From " << rank << " to "  << m_migr[0] << std::endl;
        MPI_Recv(&r_x[0], ncities, MPI_DOUBLE, m_migr[0], send02x, MPI_COMM_WORLD, &s02x);
        MPI_Send(&s_x[0], ncities, MPI_DOUBLE, m_migr[0], recv02x, MPI_COMM_WORLD);
        MPI_Recv(&r_y[0], ncities, MPI_DOUBLE, m_migr[0], send02y, MPI_COMM_WORLD, &s02y);
        MPI_Send(&s_y[0], ncities, MPI_DOUBLE, m_migr[0], recv02y, MPI_COMM_WORLD);
        MPI_Recv(&r_l[0], ncities, MPI_INTEGER, m_migr[0], send02l, MPI_COMM_WORLD, &s02l);
        MPI_Send(&s_l[0], ncities, MPI_INTEGER, m_migr[0], recv02l, MPI_COMM_WORLD);
    }


    else if (rank == m_migr[3])
    {
        std::cout <<"From " << rank << " to "  << m_migr[1] << std::endl;
        MPI_Recv(&r_x[0], ncities, MPI_DOUBLE, m_migr[1], send13x, MPI_COMM_WORLD, &s31x);
        MPI_Send(&s_x[0], ncities, MPI_DOUBLE, m_migr[1], recv13x, MPI_COMM_WORLD);
        MPI_Recv(&r_y[0], ncities, MPI_DOUBLE, m_migr[1], send13y, MPI_COMM_WORLD, &s31y);
        MPI_Send(&s_y[0], ncities, MPI_DOUBLE, m_migr[1], recv13y, MPI_COMM_WORLD);
        MPI_Recv(&r_l[0], ncities, MPI_INTEGER, m_migr[1], send13l, MPI_COMM_WORLD, &s31l);
        MPI_Send(&s_l[0], ncities, MPI_INTEGER, m_migr[1], recv13l, MPI_COMM_WORLD);


    }
/*
    if (rank == m_migr[3])
    {
        std::cout << "rank " << m_migr[3] << " receiving" << std::endl ;
        std::cout << "x = " ;
        for (int i=0; i<ncities; i++)
            std:: cout << r_x[i] << " " ;
        std::cout << std::endl ;
        std::cout << "y = " ;
        for (int i=0; i<ncities; i++)
            std:: cout << r_y[i] << " " ;
        std::cout << std::endl ;
    }
*/
    
    m_trips[0].SetPositions (r_l, r_x, r_y) ;

}
