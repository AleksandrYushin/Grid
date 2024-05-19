#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

//Распаралеливание
#include <mpi.h>

using namespace std;

// Класс нейрона
class neuron{
friend class network;
protected:
    // Параметры функции
public:
    float linear_combination(vector<float> weights, vector<float> in_neuron){
        float summ = 0;
        for (unsigned int i = 0; i < weights.size(); i++)
            summ += weights[i]*in_neuron[i];
        return summ;
    };
    virtual float operator() (float summ) = 0;
};

class neuron_sigmoid : public neuron{
protected:
    float a;
    float b;
public:
    neuron_sigmoid(float a, float b): a(a), b(b){};
    float operator() (float summ){
        return 1/pow(1+exp(-a*summ), b);
    };
};
class neuron_polynomial : public neuron{
protected:
    int n;
    float a;
    float b;
public:
    neuron_polynomial(int n, float a, float b): n(n), a(a), b(b){};
    float operator() (float summ){
        float x_lim = pow(b, 1/n) - a;
        if (summ < a)
            return 0;
        if (summ < x_lim)
            return pow(summ-a, n)/b;
        else
            return 1;
    };
};

class in_neuron : public neuron{
public:
    float operator() (float summ){
        return summ;
    };
};


// Класс сети
class network{
protected:
    vector<in_neuron*> in_neurons;
    vector<neuron*> out_neurons;
    vector<vector<neuron*>>* neurons;
    vector<vector<vector<float>>>* weights;
    vector<vector<float>>* out_weights;
public:
    network(unsigned int n, unsigned int size, unsigned int out_size){
        neurons = new vector<vector<neuron*>>;
        (*neurons).resize(n);
        weights = new vector<vector<vector<float>>>;
        (*weights).resize(n);
        for(unsigned int i = 0; i < n; i++) {
            (*neurons)[i].resize(size);
            (*weights)[i].resize(size);
            for(unsigned int j = 0; j < size; j++) {
                (*neurons)[i][j] = new neuron_sigmoid(2, 1);
                (*weights)[i][j].resize(size);
                for(unsigned int k = 0; k < size; k++) {
                    (*weights)[i][j][k] = (double)(rand())/RAND_MAX *0.02;
                };
            };
        };
        in_neurons.resize(size);
        for(unsigned int i = 0; i < size; i++) {
            in_neurons[i] = new in_neuron();
        };
        out_neurons.resize(out_size);
        out_weights = new vector<vector<float>>;
        (*out_weights).resize(out_size);
        for(unsigned int i = 0; i < out_size; i++) {
            out_neurons[i] = new neuron_sigmoid(2, 1);
            (*out_weights)[i].resize(size);
            for(unsigned int j = 0; j < size; j++) {
                (*out_weights)[i][j] = (double)(rand())/RAND_MAX *0.02;
            };
        };
    };
    ~network(){
        for(unsigned int i = 0; i < (*weights).size(); i++) {
            for(unsigned int j = 0; j < (*weights)[i].size(); j++) {
                delete (*neurons)[i][j];
            };
        };
        for(unsigned int i = 0; i < out_neurons.size(); i++) {
            delete out_neurons[i];
        };
        delete weights, out_weights, neurons;
    };
    
    double out(vector<float>* data, vector<vector<float>>* activity, vector<float>* out_activity){
        vector<float>* signal;
        vector<float>* signal_old;
        signal = new vector<float>;
        signal_old = new vector<float>;
        (*signal).resize(in_neurons.size());
        (*signal_old).resize(in_neurons.size());
        for (unsigned int i = 0; i < in_neurons.size(); i++) {
            (*signal_old)[i] = (*in_neurons[i]).operator()((*data)[i]);
        };

        (*activity).resize((*neurons).size());
        for(unsigned int i = 0; i < (*neurons).size(); i++) {
            (*activity)[i].resize((*neurons)[i].size());
            for (unsigned int j = 0; j < (*neurons)[i].size(); j++) {
                (*signal)[j] = (*neurons)[i][j]->operator()((*neurons)[i][j]->linear_combination((*weights)[i][j], *signal_old));
                (*activity)[i][j] = (*signal)[j];
            };
            (*signal_old) = *signal;
        };

        (*out_activity).resize(out_neurons.size());
        for (unsigned int i = 0; i < out_neurons.size(); i++) {
            (*signal)[i] = out_neurons[i]->operator()(out_neurons[i]->linear_combination((*out_weights)[i], *signal_old));
            (*out_activity)[i] = (*signal)[i];
        };

        float a = (*signal)[0];
        delete signal;
        delete signal_old;
        return a;
    };
    
    void get_weights(string filename){
        ofstream os(filename, ios::app);
        
        os << "weights: " << endl;
        for(unsigned int i = 0; i < (*weights).size(); i++) {
            for(unsigned int j = 0; j < (*weights)[i].size(); j++) {
                for(unsigned int k=0; k< (*weights)[i][j].size(); k++) {
                    os << (*weights)[i][j][k] << endl;
                };
            };
        };
        os << "weights from out: " << endl;
        for(unsigned int i = 0; i < (*out_weights).size(); i++) {
            for(unsigned int j = 0; j < (*out_weights)[i].size(); j++) {
                os << (*out_weights)[i][j] << endl;
            };
        };
    };
    void imput_weights(string filename){
        ifstream file(filename);
        string line;
        vector<float>* data;
        data = new vector<float>;
        float a;
        while (getline(file, line)) {
            if (line != "weights: " && line != "weights from out: ") { // если найдено кодовое слово
                istringstream iss(line);
                iss >> a;
                (*data).push_back(a);
            }
        }

        //Форматирование
        for(unsigned int i = 0; i < (*weights).size(); i++) {
            for(unsigned int j = 0; j < (*weights)[i].size(); j++) {
                for(unsigned int k=0; k< (*weights)[i][j].size(); k++) {
                    (*weights)[i][j][k] = (*data)[i*(*weights)[i].size() + j*(*weights)[i][j].size() + k];
                };
            };
        };
        int tor = (*weights).size()*(*weights)[0].size()*(*weights)[0][0].size();

        for(unsigned int i = 0; i < (*out_weights).size(); i++) {
            for(unsigned int j = 0; j < (*out_weights)[i].size(); j++) {
                (*out_weights)[i][j] = (*data)[tor +i*(*out_weights)[i].size() + j];
            };
        };

    };

    void reverse_propagation_error(float error, vector<vector<float>>* activity, vector<float>* out_activity){
        vector<float>* delta;
        delta = new vector<float>;
        (*delta).resize(in_neurons.size());
        vector<float>* delta_old;
        delta_old = new vector<float>;
        (*delta_old).resize(in_neurons.size());
        vector<float>* technic;
        technic = new vector<float>;
        (*technic).resize(in_neurons.size());

        for (unsigned int i = 0; i < out_neurons.size(); i++) {
            if (i==0){
                (*delta_old)[i] = - error*(*out_activity)[i]*(1-(*out_activity)[i]);
                for (unsigned int j = 0; j < (*neurons)[i].size(); j++){
                    (*out_weights)[i][j] -= 0.001*(*delta_old)[i]*(*activity)[(*neurons).size()-1][j];
                }
            }
            else
                (*delta_old)[i] = 0;
        };

        for(unsigned int i = (*neurons).size()-1; i > 0; i--) {
            for (unsigned int j = 0; j < (*neurons)[i].size(); j++) {
                for(unsigned int k=0; k< (*weights)[i][j].size(); k++) {
                    (*technic)[k] -= (*weights)[i][k][j];
                };
                (*delta)[j] = (*activity)[i][j]*(1-(*activity)[i][j])*((*neurons)[i][j]->linear_combination(*technic, *delta_old));
                for(unsigned int k=0; k< (*weights)[i][j].size(); k++) {
                    (*weights)[i][j][k] -= 0.001*(*delta)[j]*(*activity)[i-1][k];
                };
            };
        };

        for (unsigned int j = 0; j < (*neurons)[0].size(); j++) {
            for(unsigned int k=0; k< (*weights)[0][j].size(); k++) {
                (*technic)[k] -= (*weights)[0][k][j];
            };
            (*delta)[j] = (*activity)[0][j]*(1-(*activity)[0][j])*((*neurons)[0][j]->linear_combination(*technic, *delta_old));
            for(unsigned int k=0; k< (*weights)[0][j].size(); k++) {
                (*weights)[0][j][k] -= 0.001*(*delta)[j];
            };
        };
    };
};

vector<float> pngToVector(const string& filename) {
    int width, height, channels;
    stbi_uc* image = stbi_load(filename.c_str(), &width, &height, &channels, 0);
    
    if (!image) {
        // обработка ошибки
        return {};
    }
    
    vector<float> floatVec;
    floatVec.reserve(width * height * channels);

    for (int i = 0; i < width * height * channels; ++i) {
        floatVec.push_back(static_cast<float>(image[i]) / 255.0f);
    }

    stbi_image_free(image);
    return floatVec;
}

string intToStringWithLeadingZeros(int value, int width) {
    ostringstream oss;
    oss << std::setfill('0') << std::setw(width) << value;
    return oss.str();
}

int main(int argc, char **argv){ //Параметры распаралеливания
    network my(10, 400, 400); //49152    
    // Инициализация подсистемы MPI
    MPI_Init(&argc, &argv);

    //Обучение
    for (unsigned int k = 0; k < 1; k++){
        for (unsigned int i = 0; i < 500; i++){
           for (unsigned int j = 0; j < 10; j++){
                string name = "train_3" + to_string(j) +  "_" + intToStringWithLeadingZeros(i, 5) + ".png";
                vector<float>* data;
                data = new vector<float>;
                (*data) = pngToVector(name);

                vector<vector<float>>* activ;
                activ = new vector<vector<float>>;


                vector<float>* out_activ;
                out_activ = new vector<float>;
      
                float out = my.out(data, activ, out_activ);
        
                //Проверка вывода
                cout << endl << out << " " << float(j)/10 << endl;
                my.reverse_propagation_error(float(j)/10 - out, activ, out_activ);
                delete activ, data;
            };
        };
    };

    MPI_Finalize();
    //Вывод весов
    my.get_weights("nort.txt");
    return 0;
};


// int main(int argc, char **argv){ //Параметры распаралеливания
//     // Инициализация подсистемы MPI
//     network my(10, 400, 1);
//     my.imput_weights("nort.txt");

//     string name;
//     vector<float>* data;
//     data = new vector<float>;
//     vector<vector<float>>* activ;
//     activ = new vector<vector<float>>;

//     MPI_Init(&argc, &argv);
//     while (true)
//     {
//         cout << "NAME?" << endl;
//         cin >> name;
//         if (name == "0")
//             break;
//         (*data) = pngToVector(name);
//         float out = my.out(data, activ);
//         cout << out << endl;

//     }
//     delete activ, data;
//     MPI_Finalize();
//     return 0;
// };
