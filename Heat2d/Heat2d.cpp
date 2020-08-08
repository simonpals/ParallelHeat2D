#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <thread>
#include <fstream>
#include <string>

struct configuration {
    double BOTTOM_TEMP = 0.;
    double TOP_TEMP = 0.;
    double LEFT_TEMP = 0.;
    double RIGHT_TEMP = 0.;
    double INTERIOR_NODES = 0.;
    int ro = 0;
    double C_ro = 0.;
    double k = 0.;
    double dx = 0;
    double dy = 0;
    double dt = 0.;
    int TIME_COUNT = 0;
    int WRITE_EACH = 0;
    int THREADS_COUNT = 0;
};


int read_configurations(std::string *file_name, configuration *conf) {

    std::ifstream infile(*file_name);
    if (!infile.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return -1;
    }

    std::string BOTTOM_TEMP;
    std::string TOP_TEMP;
    std::string LEFT_TEMP;
    std::string RIGHT_TEMP;
    std::string INTERIOR_NODES;
    std::string ro;
    std::string C_ro;
    std::string k;
    std::string dx;
    std::string dy;
    std::string dt;
    std::string TIME_COUNT;
    std::string WRITE_EACH;
    std::string THREADS_COUNT;

    getline(infile, BOTTOM_TEMP);
    getline(infile, TOP_TEMP);
    getline(infile, LEFT_TEMP);
    getline(infile, RIGHT_TEMP);
    getline(infile, INTERIOR_NODES);
    getline(infile, ro);
    getline(infile, C_ro);
    getline(infile, k);
    getline(infile, dx);
    getline(infile, dy);
    getline(infile, dt);
    getline(infile, TIME_COUNT);
    getline(infile, WRITE_EACH);
    getline(infile, THREADS_COUNT);


    try {
        conf->BOTTOM_TEMP = stod(BOTTOM_TEMP);
        conf->TOP_TEMP = stod(TOP_TEMP);
        conf->LEFT_TEMP = stod(LEFT_TEMP);
        conf->RIGHT_TEMP = stod(RIGHT_TEMP);
        conf->INTERIOR_NODES = stod(INTERIOR_NODES);
        conf->ro = stoi(ro);
        conf->C_ro = stod(C_ro);
        conf->k = stod(k);
        conf->dx = stod(dx);
        conf->dy = stod(dy);
        conf->dt = stod(dt);
        conf->TIME_COUNT = stoi(TIME_COUNT);
        conf->WRITE_EACH = stoi(WRITE_EACH);
        conf->THREADS_COUNT = stoi(THREADS_COUNT);
    }
    catch (std::invalid_argument &e) {
        std::cerr << "Error, not correct configuration" << std::endl;
        return -1;
    }

    return 0;
}


void writeFileAndConsole(std::ostream &out, std::vector<std::vector<double>> &grid) {
    out << std::endl;

    for (std::vector<double> &v : grid) {
        for (double d : v) {
            out << std::setw(12) << d << " ";
        }
        out << std::endl;
    }
    out << std::endl;
}


void heat2d(configuration config, std::vector<std::vector<double>> &grid, std::vector<std::vector<double>> &grid2) {
    std::vector<std::thread> threads;
    std::vector<int> indices;

    double DIM = sqrt(config.INTERIOR_NODES) + 2;
    double alpha = config.k / (config.ro * config.C_ro);

    indices.push_back(1);

    for (int n = 0; n < config.THREADS_COUNT; n++) {
        if (n == config.THREADS_COUNT - 1) {
            indices.push_back(DIM - 1);
        }
        else {
            indices.push_back(indices[n] + DIM / config.THREADS_COUNT);
        }
    }


    for (int n = 1; n <= config.THREADS_COUNT; n++) {
        threads.emplace_back([n, &indices, &grid2, &grid, &config, &DIM, &alpha]() {
            int beg = indices[n - 1];
            int end = indices[n];
            for (int i = DIM - 2; i > 0; i--) {
                if (beg <= i && i < end) {
                    for (int j = 1; j < DIM - 1; j++) {
                        grid2[i][j] = grid[i][j] + config.dt * alpha * (((grid[i - 1][j] + -2.0 * grid[i][j]
                            + grid[i + 1][j]) / (config.dx * config.dx)) +
                            ((grid[i][j - 1] +
                                grid[i][j + 1] - 2.0 * grid[i][j]) /
                                (config.dy * config.dy)));

                    }
                }
            }
        });
    }
    for (int i = 0; i < config.THREADS_COUNT; i++) {
        threads[i].join();
    }

    grid.swap(grid2);
   //     grid = grid2;
}

int main(int argc, char *argv[]) {
    configuration config;
    std::string file_name = "configuration.txt";
    read_configurations(&file_name, &config);
    if (config.THREADS_COUNT < 1) {
        std::cerr << "Error, incorrect number of threads" << std::endl;
        return -1;
    }

    double DIM = sqrt(config.INTERIOR_NODES) + 2;
    std::vector<std::vector<double>> rows(DIM);
    std::vector<std::vector<double>> rows2(DIM);

    for (int i = 0; i < DIM; i++) {
        if (i == 0) {
            for (int j = 0; j < DIM; j++) {
                rows[i].push_back(config.TOP_TEMP);
                rows2[i].push_back(config.TOP_TEMP);
            }
        }
        else if (i > 0 && i < DIM - 1) {
            rows[i].push_back(config.LEFT_TEMP);
            rows2[i].push_back(config.LEFT_TEMP);
            for (int j = 1; j < DIM - 1; j++) {
                rows[i].push_back(0);
                rows2[i].push_back(0);
            }
            rows[i].push_back(config.RIGHT_TEMP);
            rows2[i].push_back(config.RIGHT_TEMP);
        }
        else {
            for (int j = 0; j < DIM; j++) {
                rows[i].push_back(config.BOTTOM_TEMP);
                rows2[i].push_back(config.BOTTOM_TEMP);
            }
        }
    }

    std::ofstream result("res.txt");
    std::ostream &out = (result.is_open() ? std::cout : result);

    std::cout << "=============Time: " << 0 << "===============" << std::endl;
    writeFileAndConsole(out, rows);

    int timeSteps = 0;
    do {
        heat2d(config, rows, rows2);
        timeSteps++;
        if (timeSteps % config.WRITE_EACH == 0) {
            std::cout << "=============Time: " << timeSteps * config.dt << "===============" << std::endl;
            writeFileAndConsole(out, rows);

            std::thread file_thread(writeFileAndConsole, std::ref(result), std::ref(rows));
            file_thread.join();
        }
    } while (timeSteps < config.TIME_COUNT);
    result.close();
    return 0;
}