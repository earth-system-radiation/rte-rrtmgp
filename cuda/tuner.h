#ifndef TUNER_H
#define TUNER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>


// #ifdef __CUDACC__
using Tuner_map = std::map<std::string, std::pair<dim3, dim3>>;
// #endif


class Tuner
{
    public:
        static Tuner& get()
        {
            static Tuner tuner("rte_rrtmgp_kernel_tuning.txt");
            return tuner;
        }

        ~Tuner()
        {
            std::ofstream outfile("rte_rrtmgp_kernel_tuning.txt");
            for (const auto& i : tuner_map)
            {
                const std::string del(" ");
                outfile << i.first << del << i.second.first.x << del << i.second.first.y << del << i.second.first.z
                                   << del << i.second.second.x << del << i.second.second.y << del << i.second.second.z << std::endl;
            }
            outfile.close();
        }

        static Tuner_map& get_map()
        {
            return Tuner::get().tuner_map;
        }

        std::pair<dim3, dim3> operator[](const std::string& name)
        {
            return tuner_map.at(name);
        }

    private:
        Tuner(const std::string& filename)
        {
            std::ifstream infile("rte_rrtmgp_kernel_tuning.txt");
            std::string name;
            std::array<unsigned int, 6> vals;
            while (infile >> name >> vals[0] >> vals[1] >> vals[2] >> vals[3] >> vals[4] >> vals[5])
            {
                tuner_map[name] = std::make_pair(dim3(vals[0], vals[1], vals[2]), dim3(vals[3], vals[4], vals[5]));
            }
            infile.close();
        }

        Tuner(const Tuner&) = delete;
        Tuner& operator=(const Tuner&) = delete;

        Tuner_map tuner_map;
};


inline bool file_exists(const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}


////  RUNTIME TUNER ////
template<class Func, class... Args>
std::tuple<dim3, dim3> tune_kernel(
        const std::string& kernel_name,
        dim3 problem_size,
        const std::vector<unsigned int>& ib, const std::vector<unsigned int>& jb, const std::vector<unsigned int>& kb,
        Func&& f, Args&&... args)
{
    std::cout << "Tuning " << kernel_name << ": ";

    std::ofstream tuner_output;

    std::string file_name = kernel_name;
    if (file_exists(file_name + ".txt"))
    {
        for (int i=2;; ++i)
        {
            if (!file_exists(file_name + "_" + std::to_string(i) + ".txt"))
            {
                file_name += "_" + std::to_string(i);
                break;
            }
        }
    }

    tuner_output.open(file_name + ".txt", std::ios::out | std::ios::app);
    tuner_output
        << std::setw(10) << "block_x"
        << std::setw(10) << "block_y"
        << std::setw(10) << "block_z"
        << std::setw(16) << "time (ms)" << "\n";

    float fastest = std::numeric_limits<float>::max();
    dim3 fastest_block{1, 1, 1};
    dim3 fastest_grid{problem_size};

    for (const unsigned int k : kb)
        for (const int unsigned j : jb)
            for (const unsigned int i : ib)
            {
                dim3 block{i, j, k};
                dim3 grid{
                    problem_size.x/block.x + (problem_size.x%block.x > 0),
                    problem_size.y/block.y + (problem_size.y%block.y > 0),
                    problem_size.z/block.z + (problem_size.z%block.z > 0)};

                constexpr int n_samples = 8;

                // Warmup...
                for (int n=0; n<n_samples; ++n)
                    f<<<grid, block>>>(args...);

                cudaDeviceSynchronize();
                cudaEvent_t start;
                cudaEvent_t stop;
                cudaEventCreate(&start);
                cudaEventCreate(&stop);

                cudaEventRecord(start, 0);
                for (int n=0; n<n_samples; ++n)
                    f<<<grid, block>>>(args...);
                cudaEventRecord(stop, 0);

                cudaEventSynchronize(stop);
                float duration = 0.f;
                cudaEventElapsedTime(&duration, start, stop);

                cudaEventDestroy(start);
                cudaEventDestroy(stop);

                // Check whether kernel has succeeded.
                cudaError err = cudaGetLastError();
                if (err != cudaSuccess)
                {
                    tuner_output
                        << std::setw(10) << i
                        << std::setw(10) << j
                        << std::setw(10) << k
                        << std::setw(16) << "NaN\n";
                }
                else
                {
                    if (duration < fastest)
                    {
                        fastest = duration;
                        fastest_grid = grid;
                        fastest_block = block;
                    }

                    tuner_output
                        << std::setw(10) << i
                        << std::setw(10) << j
                        << std::setw(10) << k
                        << std::setw(16) << duration/n_samples << "\n";
                }
            }

    tuner_output.close();

    std::cout << "("
        << fastest_block.x << ", "
        << fastest_block.y << ", "
        << fastest_block.z << ")" << std::endl;

    return {fastest_grid, fastest_block};
}


////  COMPILE TIME TUNER ////
// Tuning functions.
template<class Func, unsigned int I, unsigned int J, unsigned int K, class... Args>
void tune_ijk(
        dim3 problem_size, dim3& fastest_grid, dim3& fastest_block, float& fastest,
        std::ofstream& tuner_output,
        Args... args)
{
    dim3 block{I, J, K};
    dim3 grid{
        problem_size.x/block.x + (problem_size.x%block.x > 0),
        problem_size.y/block.y + (problem_size.y%block.y > 0),
        problem_size.z/block.z + (problem_size.z%block.z > 0)};

    // Warmup...
    constexpr int n_samples = 8;

    for (int i=0; i<n_samples; ++i)
        Func::template launch<I, J, K>(grid, block, args...);

    cudaDeviceSynchronize();
    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);
    for (int i=0; i<n_samples; ++i)
        Func::template launch<I, J, K>(grid, block, args...);
    cudaEventRecord(stop, 0);

    cudaEventSynchronize(stop);
    float duration = 0.f;
    cudaEventElapsedTime(&duration, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // Check whether kernel has succeeded.
    cudaError err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        tuner_output
            << std::setw(10) << I
            << std::setw(10) << J
            << std::setw(10) << K
            << std::setw(16) << "NaN\n";
    }
    else
    {
        if (duration < fastest)
        {
            fastest = duration;
            fastest_grid = grid;
            fastest_block = block;
        }

        tuner_output
            << std::setw(10) << I
            << std::setw(10) << J
            << std::setw(10) << K
            << std::setw(16) << duration/n_samples << "\n";
    }
}


template<class Func, unsigned int I, unsigned int J, unsigned int... Ks, class... Args>
void tune_ij(
        std::integer_sequence<unsigned int, Ks...> ks,
        dim3 problem_size, dim3& fastest_grid, dim3& fastest_block, float& fastest,
        std::ofstream& tuner_output,
        Args... args)
{
    (tune_ijk<Func, I, J, Ks>(problem_size, fastest_grid, fastest_block, fastest, tuner_output, args...), ...);
}


template<class Func, unsigned int I, unsigned int... Js, unsigned int... Ks, class... Args>
void tune_i(
        std::integer_sequence<unsigned int, Js...> js,
        std::integer_sequence<unsigned int, Ks...> ks,
        dim3 problem_size, dim3& fastest_grid, dim3& fastest_block, float& fastest,
        std::ofstream& tuner_output,
        Args... args)
{
    (tune_ij<Func, I, Js>(ks, problem_size, fastest_grid, fastest_block, fastest, tuner_output, args...), ...);
}


template<class Func, class... Args, unsigned int... Is, unsigned int... Js, unsigned int... Ks>
std::tuple<dim3, dim3> tune_kernel_compile_time(
        const std::string& kernel_name,
        dim3 problem_size,
        std::integer_sequence<unsigned int, Is...> is,
        std::integer_sequence<unsigned int, Js...> js,
        std::integer_sequence<unsigned int, Ks...> ks,
        Args&&... args)
{
    std::cout << "Tuning " << kernel_name << ": ";

    std::ofstream tuner_output;

    std::string file_name = kernel_name;
    if (file_exists(file_name + ".txt"))
    {
        for (int i=2;; ++i)
        {
            if (!file_exists(file_name + "_" + std::to_string(i) + ".txt"))
            {
                file_name += "_" + std::to_string(i);
                break;
            }
        }
    }

    tuner_output.open(file_name + ".txt", std::ios::out);
    tuner_output
        << std::setw(10) << "block_x"
        << std::setw(10) << "block_y"
        << std::setw(10) << "block_z"
        << std::setw(16) << "time (ms)" << "\n";

    float fastest = std::numeric_limits<float>::max();
    dim3 fastest_block{1, 1, 1};
    dim3 fastest_grid{problem_size};

    (tune_i<Func, Is>(js, ks, problem_size, fastest_grid, fastest_block, fastest, tuner_output, args...), ...);

    std::cout << "("
        << fastest_block.x << ", "
        << fastest_block.y << ", "
        << fastest_block.z << ")" << std::endl;

    return {fastest_grid, fastest_block};
}


// Running functions.
template<class Func, unsigned int I, unsigned int J, unsigned int K, class... Args>
void run_ijk(
        dim3 grid, dim3 block,
        Args... args)
{
    if ( (block.x == I) && (block.y == J) && (block.z == K) )
        Func::template launch<I, J, K>(grid, block, args...);
}


template<class Func, unsigned int I, unsigned int J, unsigned int... Ks, class... Args>
void run_ij(
        std::integer_sequence<unsigned int, Ks...> ks,
        dim3 grid, dim3 block,
        Args... args)
{
    (run_ijk<Func, I, J, Ks>(grid, block, args...), ...);
}


template<class Func, unsigned int I, unsigned int... Js, unsigned int... Ks, class... Args>
void run_i(
        std::integer_sequence<unsigned int, Js...> js,
        std::integer_sequence<unsigned int, Ks...> ks,
        dim3 grid, dim3 block,
        Args... args)
{
    (run_ij<Func, I, Js>(ks, grid, block, args...), ...);
}


template<class Func, class... Args, unsigned int... Is, unsigned int... Js, unsigned int... Ks>
void run_kernel_compile_time(
        std::integer_sequence<unsigned int, Is...> is,
        std::integer_sequence<unsigned int, Js...> js,
        std::integer_sequence<unsigned int, Ks...> ks,
        dim3 grid, dim3 block,
        Args&&... args)
{
    (run_i<Func, Is>(js, ks, grid, block, args...), ...);
}
#endif
