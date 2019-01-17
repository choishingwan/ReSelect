#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <getopt.h>
#include <unordered_set>
#include <unordered_map>
#include <assert.h>
#include <sstream>
#include <math.h>
#include "gzstream.h"


template <class T>
class vec2d
{
public:
    vec2d() {}
    vec2d(size_t row, size_t col, T def)
    {
        if (row == 0 || col == 0) {
            throw std::invalid_argument("Dimension of 2d vector must be >0");
        }
        m_storage.resize(row * col, def);
        m_row = row;
        m_col = col;
    }
    vec2d(size_t row, size_t col)
    {
        if (row == 0 || col == 0) {
            throw std::invalid_argument("Dimension of 2d vector must be >0");
        }
        m_storage.resize(row * col);
        m_row = row;
        m_col = col;
    }
    T operator()(size_t row, size_t col) const
    {
        if (row > m_row || col > m_col){
            std::cerr << "In const" << std::endl;
            std::cerr << "Dimension: " << m_row << "\t" << m_col << std::endl;
            std::cerr << "Requested: " << row << "\t" << col << std::endl;
            throw std::out_of_range("2d vector out of range!");
        }
        return m_storage[row * m_col + col];
    }
    T& operator()(size_t row, size_t col)
    {
        if (row > m_row || col > m_col){
            std::cerr << "Dimension: " << m_row << "\t" << m_col << std::endl;
            std::cerr << "Requested: " << row << "\t" << col << std::endl;
            throw std::out_of_range("2d vector out of range!");
        }
        return m_storage[row * m_col + col];
    }
    void clear() { m_storage.clear(); }
    size_t rows() const { return m_row; }
    size_t cols() const { return m_col; }

private:
    size_t m_row = 0;
    size_t m_col = 0;
    std::vector<T> m_storage;
};



void usage()
{

    std::cerr << "Relatedness Selector\n";
    std::cerr << "Usage: ReSelect [options]\n";
    std::cerr << "Options:\n";
    std::cerr << "    --input   | -i    Relatedness matrix file\n";
    std::cerr << "    --id      | -d    Relatedness ID file\n";
    std::cerr << "    --base    | -b    Base Sample IDs\n";
    std::cerr << "    --target  | -t    Target Sample IDs\n";
    std::cerr << "    --out     | -o    Output Prefix\n";
    std::cerr << "    --help    | -h    Display this help message\n";
}

std::vector<std::string> split(const std::string& seq,
                               const std::string& separators=" \t")
{
    std::size_t prev = 0, pos;
    std::vector<std::string> result;
    while ((pos = seq.find_first_of(separators, prev)) != std::string::npos)
    {
        if (pos > prev) result.emplace_back(seq.substr(prev, pos - prev));
        prev = pos + 1;
    }
    if (prev < seq.length())
        result.emplace_back(seq.substr(prev, std::string::npos));
    return result;
}
void ltrim(std::string& s)
{
    s.erase(s.begin(),
            std::find_if(s.begin(), s.end(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))));
}
// trim from end (in place)
void rtrim(std::string& s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace)))
                .base(),
            s.end());
}
// trim from both ends (in place)
void trim(std::string& s)
{
    ltrim(s);
    rtrim(s);
}

// function from John D.Cook
// https://www.johndcook.com/blog/standard_deviation/
class RunningStat
{
public:
    RunningStat() {}
    void clear()
    {
        n = 0;
        M1 = M2 = M3 = M4 = 0.0;
    }
    void push(double x)
    {
        double delta, delta_n, delta_n2, term1;

        size_t n1 = n;
        n++;
        delta = x - M1;
        assert(n > 0);
        delta_n = delta / n;
        delta_n2 = delta_n * delta_n;
        term1 = delta * delta_n * n1;
        M1 += delta_n;
        M4 += term1 * delta_n2 * (n * n - 3 * n + 3) + 6 * delta_n2 * M2
              - 4 * delta_n * M3;
        M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
        M2 += term1;
    }
    size_t get_n() const { return n; }

    double mean() const { return M1; }

    double var() const { return M2 / ((double) n - 1.0); }

    double sd() const { return sqrt(var()); }

private:
    size_t n = 0;
    double M1 = 0, M2 = 0, M3 = 0, M4 = 0;
};


template <typename T>
inline T convert(const std::string& str)
{
    std::istringstream iss(str);
    T obj;
    iss >> obj;

    if (!iss.eof() || iss.fail()) {
        throw std::runtime_error("Unable to convert the input");
    }
    return obj;
}


int main(int argc, char *argv[])
{
    if(argc < 2){
        usage();
        return -1;
    }
    static const char* optString = "i:d:b:t:o:h?";
    static const struct option longOpts[] = {
    {"input", required_argument, nullptr, 'i'},
    {"id", required_argument, nullptr, 'd'},
    {"base", required_argument, nullptr, 'b'},
    {"target", required_argument, nullptr, 't'},
    {"out", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

    int longIndex = 0;
    int opt = 0;
    std::string command = "";
    opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    std::string error_message = "";

    std::string base_file_name;
    std::string id_file_name;
    std::string matrix_file_name;
    std::string out_name;
    std::string target_file_name;

    while (opt != -1) {
        switch (opt)
        {
        case 'i': matrix_file_name = optarg; break;
        case 'o': out_name = optarg; break;
        case 'b': base_file_name = optarg; break;
        case 'd': id_file_name = optarg; break;
        case 't': target_file_name = optarg; break;
        case 'h':
        case '?':
            usage();
            exit(0);
        default:
            throw "Undefined operator, please use --help for more information!";
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    bool error = false;
    if(matrix_file_name.empty()){
        error =true;
        std::cerr << "Error: You must provide the relationship matrix file" << std::endl;
    }
    if(id_file_name.empty()){
        error = true;
        std::cerr << "Error: You must provide the relationship matrix ID file" << std::endl;
    }
    if(base_file_name.empty()){
        error = true;
        std::cerr << "Error: You must provide a file containing sample IDs in base" << std::endl;
    }
    if(target_file_name.empty() && !error){
        target_file_name = base_file_name;
        std::cerr << "Warning: File containing sample IDs in target not provided. Will use file for base" << std::endl;
    }
    if(error){
        std::cerr << "Please check you have the correct inputs" << std::endl;
        return -1;
    }
    // 1 means it is row, 2 means it is col, 3 means it is both, 0 means it is none
    std::unordered_set<std::string> base_ids;
    std::unordered_set<std::string> target_ids;

    std::unordered_map<std::string, int> sample_id;
    std::ifstream base_file, target_file, matrix_file, id_file;
    GZSTREAM_NAMESPACE::igzstream matrix_gz_file;
    bool gz_input = false;
    if (matrix_file_name.substr(matrix_file_name.find_last_of(".") + 1).compare("gz") == 0) {
        matrix_gz_file.open(matrix_file_name.c_str());
        if (!matrix_gz_file.good()) {
            std::string error_message = "Error: Cannot open GRM file: "
                                        + matrix_file_name + " (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
        gz_input = true;
    }
    id_file.open(id_file_name.c_str());
    if(!id_file.is_open()){
        error =true;
        std::cerr << "Error: Cannot open relationship matrix ID file - " << id_file_name << std::endl;
    }
    base_file.open(base_file_name.c_str());
    if(!base_file.is_open()){
        error = true;
        std::cerr << "Error: Cannot open base ID file - " << base_file_name << std::endl;
    }
    if(base_file_name != target_file_name){
        target_file.open(target_file_name.c_str());
        if(!target_file.is_open()){
            error = true;
            std::cerr << "Error: Cannot open target ID file - " << target_file_name << std::endl;
        }
    }
    matrix_file.open(matrix_file_name.c_str());
    if(!matrix_file.is_open()){
        error = true;
        std::cerr << "Error: Cannot open relationship matrix file - " << matrix_file_name << std::endl;
    }

    std::ofstream output_matrix;
    output_matrix.open(std::string(out_name+".matrix").c_str());
    if(!output_matrix.is_open()){
        error = true;
        std::cerr << "Cannot open output file to write - " << std::string(out_name+".matrix") << std::endl;
    }
    std::ofstream output_avg;
    output_avg.open(std::string(out_name+".avg").c_str());
    if(!output_avg.is_open()){
        error = true;
        std::cerr << "Cannot open output file to write - " << std::string(out_name+".avg") << std::endl;
    }
    std::ofstream output_row;
    output_row.open(std::string(out_name+".row").c_str());
    if(!output_row.is_open()){
        error = true;
        std::cerr << "Cannot open output file to write - " << std::string(out_name+".row") << std::endl;
    }
    std::ofstream output_col;
    output_col.open(std::string(out_name+".col").c_str());
    if(!output_col.is_open()){
        error = true;
        std::cerr << "Cannot open output file to write - " << std::string(out_name+".col") << std::endl;
    }
    if(error){
        std::cerr << "Please check you have the correct file(s)" << std::endl;
        return -1;
    }
    std::string line;
    std::vector<std::string> token;
    std::string id;
    std::cerr << "Reading Base ID file" << std::endl;
    while(std::getline(base_file, line)){
        trim(line);
        if(line.empty()) continue;
        token = split(line);
        id = token[0] + "_" + token[1];
        base_ids.insert(id);
    }
    base_file.close();
    std::cerr << base_ids.size() << " samples in Base ID file" << std::endl;
    if(base_file_name != target_file_name){
        std::cerr << "Reading Target ID file" << std::endl;
        assert(target_file.is_open());
        while(std::getline(target_file, line)){
            trim(line);
            if(line.empty()) continue;
            token = split(line);
            id = token[0] + "_" + token[1];
            target_ids.insert(id);
        }
        target_file.close();
    }else{
        target_ids = base_ids;
    }
    std::cerr << target_ids.size() << " samples in Target ID file" << std::endl;
    // now start generating the result file at the same time as reading the file

    std::vector<bool> in_base, in_target;
    std::vector<int> base_idx, target_idx;
    int base_idx_it = 0, target_idx_it = 0;
    size_t num_base=0, num_target = 0;
    bool base, target;
    while(getline(id_file, line)){
        trim(line);
        if(line.empty()) continue;
        base= false;
        target = false;
        token = split(line);
        id = token[0]+"_"+token[1];
        base = (base_ids.find(id)!=base_ids.end());
        target = (target_ids.find(id)!=target_ids.end());
        in_base.push_back(base);
        in_target.push_back(target);
        if(base){
            base_idx.push_back(base_idx_it++);
            output_row << token[0] << "\t" << token[1] << std::endl;
            ++num_base;
        }else{
            base_idx.push_back(-1);
        }
        if(target){
            target_idx.push_back(target_idx_it++);
            output_col << token[0] << "\t" << token[1] << std::endl;
            ++num_target;
        }else{
            target_idx.push_back(-1);
        }
    }
    std::cerr << num_base << " Base ID found in relationship matrix" << std::endl;
    std::cerr << num_target << " Target ID found in relationship matrix" << std::endl;
    vec2d<double> rel_matrix(num_base, num_target);
    base_ids.clear();
    target_ids.clear();
    output_row.close();
    output_col.close();
    id_file.close();
    size_t cur_id = 0;
    while((!gz_input && std::getline(matrix_file, line))
          || (gz_input && std::getline(matrix_gz_file, line))){
        trim(line);
        if(line.empty()) continue;
        if(in_base[cur_id] || in_target[cur_id]){
            // only do the printing and calculation if we need this row
            token = split(line);
            for(size_t i = 0; i < token.size(); ++i){
                double ibd = convert<double>(token[i]);
                if(i == cur_id) ibd = 1.0;
                if(in_base[cur_id] && in_target[i]){
                    rel_matrix(static_cast<size_t>(base_idx[cur_id]),
                               static_cast<size_t>(target_idx[i])) = ibd;
                }
                if(in_base[i] && in_target[cur_id]){
                    rel_matrix(static_cast<size_t>(base_idx[i]),
                               static_cast<size_t>(target_idx[cur_id])) = ibd;
                }
            }
        }
        ++cur_id;
    }
    if(!gz_input) matrix_file.close();
    else matrix_gz_file.close();

    RunningStat mean_related;
    for(size_t r = 0; r < rel_matrix.rows(); ++r){
        mean_related.push(rel_matrix(r, 0));
        output_matrix<< rel_matrix(r, 0);
        for(size_t c = 1; c < rel_matrix.cols(); ++c){
            output_matrix << "\t" << rel_matrix(r, c);
            mean_related.push(rel_matrix(r, c));
        }
        output_matrix << std::endl;
    }
    output_matrix.close();
    output_avg << "Mean\tSD\tVar\tN" << std::endl;
    output_avg << mean_related.mean() << "\t" << mean_related.sd() << "\t" << mean_related.var() << "\t" << mean_related.get_n() << std::endl;

    output_avg.close();
    std::cerr << "Mean Relatedness is: " <<mean_related.mean() << std::endl;
    return 0;
}
