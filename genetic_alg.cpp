#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "CImg.h"
#include <iomanip>
using namespace std;
using namespace cimg_library;

// вар 14: I(А)+II(A,C)+III(A,B)+IV(A,B,I)+V(A,F)+VI(B)+VII(А)
// ф-ция: f(x)=3x^3-2x+5,x из [1..10]

struct Const {
    inline static int xKol=10;    // кол-во допустимых значений x: 1..10
    inline static int dlHro=4;     // длина хромосомы
    inline static int razPop=10;   // размер популяции
    inline static double verKro =0.7;
    inline static double verMut=0.2;
    inline static int kolPok= 50;
    inline static int parSel=1;    // кол-во пар для селекции за одну итерацию
    inline static int left =1;      // левая граница x
};

struct Osob {
    vector<short> hrom=vector<short>(Const::dlHro,0);
    int znach= 0;
};

int Func(int x) {
    return 3*x *x* x - 2*x +5;
}

int decode(const vector<short>& hrom) {
    int a= 0;
    for (int k=0; k<Const::dlHro; k++) a= a *2 +hrom[k];
    return a+ 1; 
}

vector<short> encode(int x) {
    vector<short> h(Const::dlHro,0);
    int tmp=x-1;
    for (int j=Const::dlHro -1; j >=0; j--,tmp/= 2) h[j]=tmp%2;
    return h;
}

bool isValid(const vector<short>& hrom)
{
    int a=0;
    for (int k=0; k<Const::dlHro; k++) a= a*2 +hrom[k];
    return (a<Const::xKol);
}


void grfGr(const vector<int>& num){
    if (num.empty()) return;

    long long mnX=1;
    long long mxX=num.size();

    double mnY=0;
    double mxY=*max_element(num.begin(),num.end());
    mxY*= 1.1;

    CImg<unsigned char> img(1200, 800,1,3,255);
    const unsigned char black[]={0, 0,0};
    const unsigned char blue[]={0,0, 255};
    const unsigned char red[]={255,0, 0};
    const unsigned char gray[]={200, 200, 200};

    int otLev = 150,otNiz = 100,otVrh=60,otPrv= 50;
    int shGrf=img.width() -otLev-otPrv;
    int vsGrf=img.height()- otNiz -otVrh;

    auto trX=[&](long long x) -> int {
        double dolja=static_cast<double>(x - mnX) / (mxX -mnX);
        return otLev+static_cast<int>(dolja * shGrf);
    };

    auto trY=[&](double y) -> int {
        return img.height()-otNiz-static_cast<int>((y-mnY) * vsGrf / (mxY-mnY));
    };

    img.draw_line(otLev,img.height()- otNiz,img.width()-otPrv,img.height()-otNiz,black);
    img.draw_line(otLev,img.height() -otNiz,otLev,otVrh,black);
    img.draw_text(img.width() / 2 - 80,img.height()-40,"Generation num",black, 0, 1, 18);

    int delX=14;
    for (int i=0; i <= delX; i++) {
        long long znX=mnX + i * (mxX-mnX) / delX;
        int pikX=trX(znX);
        img.draw_line(pikX,img.height()-otNiz,pikX,img.height()-otNiz +5,black);
        stringstream ss; ss<<znX;
        img.draw_text(pikX - 40,img.height()-otNiz + 15,ss.str().c_str(),black, 0 ,1 ,12);
        img.draw_line(pikX,img.height()-otNiz,pikX,otVrh,gray);
   }

    int delY=10;
    for (int i=0; i <= delY; i++) {
        double znY=mnY+i*(mxY-mnY) / delY;
        int pikY=trY(znY);
        img.draw_line(otLev,pikY,otLev - 5,pikY,black);
        stringstream ss; ss<<static_cast<int>(znY);
        img.draw_text(otLev-145,pikY-8,ss.str().c_str(),black, 0, 1, 12);
        img.draw_line(otLev,pikY,img.width()-otPrv,pikY,gray);
   }

    for (size_t i=0; i<num.size(); i++) {
        int x=trX(i+1);
        int y=trY(num[i]);
        if (i==num.size()-1) {
            img.draw_circle(x, y, 7,red);
            img.draw_circle(x,y, 7, red, 1.0f, ~0U);
       }
        else {
            img.draw_circle(x,y,5 ,blue);
            img.draw_circle(x,y,5, blue ,1.0f ,~0U);
       }
   }

    img.draw_text(img.width() / 2- 150,10,"The value of the descendant function",black , 0, 1 ,20);
    img.save("genetic_algorithm_function_values.bmp");
    cout<<"график сохранён в genetic_algorithm_function_values.bmp"<<endl;
}


void ustPar()
{
    int vybor;
    string vvod;

    cout<<"задать размер популяции? (1-да,2-нет)\n- ";
    cin >> vybor; cin.ignore( 10000,'\n');
    if (vybor == 1) {
        do {
            cout<<"Размер популяции [10;100]:\n- ";
            getline(cin ,vvod);
            Const::razPop= stoi( vvod );
            if (Const::razPop<10 || Const::razPop > 100)
                cout<<"ошибка: значение должно быть в [10;100]\n";
       } while (Const::razPop<10|| Const::razPop > 100);
   }
    else cout<<"значение по умолчанию: "<<Const::razPop<<"\n";

    cout<<"\nзадать вер-ть кроссинговера? (1-да,2-нет)\n- ";
    cin >> vybor; cin.ignore(10000,'\n');
    if (vybor==1) {
        do {
            cout<<"вер-ть кроссинговера [0.1 ; 1.0]:\n- ";
            getline(cin,vvod);
            Const::verKro=stod(vvod);
            if (Const::verKro<0.1 || Const::verKro > 1.0)
                cout<<"ошибка: значение должно быть в [0.1;1.0]\n";
       } while (Const::verKro<0.1 || Const::verKro > 1.0);
   }
    else cout<<"значение по умолчанию: "<<Const::verKro<<"\n";

    cout<<"\nзадать вер-ть мутации? (1-да,2-нет)\n- ";
    cin >> vybor; 
    cin.ignore( 10000,'\n');
    if (vybor==1) {
        do {
            cout<<"вер-ть мутации [0.1;1.0]:\n- ";
            getline(cin, vvod);
            Const::verMut=stod( vvod);
            if (Const::verMut<0.1 || Const::verMut > 1.0)
                cout<<"ошибка: значение должно быть в [0.1;1.0]\n";
       } while (Const::verMut<0.1 || Const::verMut > 1.0);
   }
    else cout<<"значение по умолчанию: "<<Const::verMut<<"\n";

    cout<<"\nзадать число поколений? (1-да,2-нет)\n- ";
    cin >> vybor; 
    cin.ignore( 10000,'\n');
    if (vybor == 1) {
        do {
            cout<<"число поколений [10;1000]:\n- ";
            getline(cin ,vvod);
            Const::kolPok =stoi(vvod);
            if (Const::kolPok<10 || Const::kolPok > 1000)
                cout<<"ошибка: значение должно быть в [10;1000]\n";
       } while (Const::kolPok<10 || Const::kolPok > 1000);
   }
    else cout<<"значение по умолчанию: "<<Const::kolPok<<"\n";
}

void genOdj(vector<Osob>& tekPop)
{
    int n=Const::razPop;
    for (int i= 0; i<n; i++) {
        int x=1+(i%Const::xKol); 
        tekPop[i].znach = Func(x);
        tekPop[i].hrom = encode(x);
   }
}

void genFok(vector<Osob>& tekPop)
{
    int center=1+rand()%Const::xKol;
    for (int i=0; i<Const::razPop; i++) {
        int offset=rand()% 3 - 1;
        int x=center+offset;
        x=max(1,min(Const::xKol,x));
        tekPop[i].znach =Func(x);
        tekPop[i].hrom= encode(x);
   }
}

vector<pair<int,int>> selSlu()
{
    vector<pair<int,int>> parVyb(Const::parSel,make_pair(-1,-1));
    for (int i=0; i<Const::parSel; i++) {
        int per = rand()%Const::razPop;
        int vtor=rand()%Const::razPop;
        while (vtor==per) vtor=rand()%Const::razPop;
        parVyb[i]={per,vtor};
   }
    return parVyb;
}

vector<pair<int,int>> selShk(const vector<Osob>& tekPop)
{
    double sumF=0.0;
    for (const auto& osob : tekPop) sumF += osob.znach;

    if (sumF <= 0.0) return selSlu();


    vector<double> cdf(Const::razPop,0.0);
    double cumul = 0.0;
    for (int i=0; i<Const::razPop; i++) {
        cumul += (double)tekPop[i].znach / sumF;
        cdf[i]=cumul;
    }
    cdf[Const::razPop-1] = 1.0;

    vector<pair<int,int>> parVyb(Const::parSel,{-1,-1});
    for (int i=0; i<Const::parSel; i++) {
        auto spinWheel=[&]() -> int {
            double r=(double)(rand()%10000) / 10000.0;
            for (int j=0; j<Const::razPop; j++)
                if (r <= cdf[j]) return j;
            return Const::razPop-1;
           };

        int per=spinWheel();
        int vtor=spinWheel();
        int tries = 0;
        while (vtor == per && tries <100) {
            vtor =spinWheel();
            tries++;
       }
        if (vtor==per) {
            vtor=(per+1) % Const::razPop;
       }
        parVyb[i] = {per,vtor};
   }
    return parVyb;
}

void kroOdn(const vector<pair<int,int>>& rodit,const vector<Osob>& tekPop,vector<Osob>& novPop)
{
    for (const auto& par : rodit) {
        double prob=(double)(rand() % 10000) / 10000.0;
        if (prob<Const::verKro) {
            int tochRz=1+rand() % (Const::dlHro-1);

            const auto& p1 = tekPop[par.first].hrom;
            const auto& p2 = tekPop[par.second].hrom;

            vector<short> pot1(p1.begin(),p1.begin()+tochRz);
            pot1.insert(pot1.end(),p2.begin()+tochRz,p2.end());

            vector<short> pot2(p2.begin(),p2.begin()+tochRz);
            pot2.insert(pot2.end(),p1.begin()+tochRz,p1.end());

            if (isValid(pot1)) {
                int x = decode(pot1);
                novPop.push_back({pot1,Func(x)});
           }
            if (isValid(pot2)) {
                int x = decode(pot2);
                novPop.push_back({pot2,Func(x)});
           }
       }
   }
}

void kroDvu(const vector<pair<int,int>>& rodit,const vector<Osob>& tekPop,vector<Osob>& novPop)
{
    for (const auto& par : rodit) {
        double prob=(double)(rand()%10000) / 10000.0;
        if (prob<Const::verKro) {
            int tch1 = 1+ rand()%(Const::dlHro-1);
            int tch2 = 1 + rand()%(Const::dlHro-1);
            while (tch1==tch2) tch2 =1+rand() % (Const::dlHro-1);
            if (tch1 > tch2) swap(tch1,tch2);

            const auto& p1 = tekPop[par.first].hrom;
            const auto& p2 = tekPop[par.second].hrom;

            vector<short> pot1(p1.begin(),p1.begin()+tch1);
            pot1.insert(pot1.end(),p2.begin()+tch1,p2.begin()+tch2);
            pot1.insert(pot1.end(),p1.begin()+tch2,p1.end());

            vector<short> pot2(p2.begin(),p2.begin()+tch1);
            pot2.insert(pot2.end(),p1.begin()+tch1,p1.begin()+tch2);
            pot2.insert(pot2.end(),p2.begin()+tch2,p2.end());

            if (isValid(pot1)) {
                int x =decode(pot1);
                novPop.push_back({pot1,Func(x)});
           }
            if (isValid(pot2)) {
                int x =decode(pot2);
                novPop.push_back({pot2,Func(x)});
           }
       }
   }
}

void kroCik(const vector<pair<int,int>>& rodit,const vector<Osob>& tekPop,vector<Osob>& novPop)
{
    for (const auto& par : rodit) {
        double prob=(double)(rand()%10000) / 10000.0;
        if (prob<Const::verKro) {
            const auto& p1 =tekPop[par.first].hrom;
            const auto& p2 =tekPop[par.second].hrom;
            int L=Const::dlHro;

            vector<bool> inCycle(L,false);
            int cur=0;
            for (int step=0; step<L; step++) {
                inCycle[cur] = true;
                short val = p1[cur];
                int next =- 1;
                for (int j= 0; j<L; j++) {
                    if ( !inCycle[j] && p2[j]==val) {
                        next =j;
                        break;
                   }
               }
                if (next == -1) break;
                cur = next;
           }

            vector<short> pot1(L),pot2(L);
            for (int j=0; j<L; j++) {
                if (inCycle[j]) {
                    pot1[j] = p1[j];
                    pot2[j] = p2[j];
               }
                else {
                    pot1[j] = p2[j];
                    pot2[j] = p1[j];
               }
           }

            if (isValid( pot1 )) {
                int x=decode(pot1);
                novPop.push_back({pot1,Func(x)});
           }
            if (isValid( pot2)) {
                int x= decode(pot2);
                novPop.push_back({pot2,Func(x)});
           }
       }
   }
}

Osob mutPro(Osob osob)
{
    double prob = (double) (rand()%10000) / 10000.0;
    if (prob<Const::verMut) {
        int ind = rand() %Const::dlHro;
        osob.hrom[ind] ^= 1;
        if (!isValid(osob.hrom)) {
            osob.hrom[ind] ^= 1;
       }
        else {
            osob.znach=Func(decode(osob.hrom));
       }
   }
    return osob;
}

Osob mutInv(Osob osob)
{
    double prob=(double)(rand()%10000) / 10000.0;
    if (prob<Const::verMut) {
        int pos1 = rand() % Const::dlHro;
        int pos2 = rand() % Const::dlHro;
        if (pos1 > pos2) swap(pos1,pos2);
        vector<short> backup=osob.hrom;
        reverse(osob.hrom.begin() + pos1,osob.hrom.begin()+pos2+1);
        if (!isValid(osob.hrom)) {
            osob.hrom=backup;
       }
        else {
            osob.znach=Func(decode(osob.hrom));
       }
   }
    return osob;
}

void applyMutation(vector<Osob>& novPop,int mut_choice)
{
    for (auto& osob : novPop) {
        if (mut_choice == 1) osob=mutPro(osob);
        else if (mut_choice == 2) osob=mutInv(osob);
        else {
            if (rand() % 2 == 0) osob=mutPro(osob);
            else osob=mutInv(osob);
       }
   }
}

Osob selElt(const vector<Osob>& novPop)
{
    if (novPop.empty()) return Osob{};
    int indNg = 0 ;
    for (int i=1; i<(int)novPop.size() ; ++i)
        if (novPop[i].znach > novPop[indNg].znach)
            indNg=i ;
    return novPop[indNg];
}

int main()
{
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
    srand((unsigned)time(0)) ;

    ustPar();

    vector<Osob> tekPop(Const::razPop) ;
    vector<Osob> novPop;

    do {
        int ind =-1;
        cout<<"стратегия начальной популяции: 1-\"одеяло\" (II.A),2-фокусировка (II.C)\n- ";
        cin >> ind;
        cin.ignore(10000,'\n');

        if (ind == 1) {
            genOdj(tekPop); 
            break; 
       }

        if (ind == 2) {
            genFok(tekPop); 
            break; 
       }

        cout<<"неверный ввод."<<endl;
   } while (true);

    int cross_choice=1;
    do {
        cout<<"оператор кроссинговера: 1-Одноточечный (IV.A),2-Двухточечный (IV.B),3-Циклический (IV.I)\n- ";
        cin >> cross_choice;
        cin.ignore( 10000,'\n');
        if (cross_choice >= 1 && cross_choice <= 3) break;
        cout<<"неверный ввод."<<endl;
   } while (true);

    int select_choice=1;
    do {
        cout<<"селекция родителей: 1-Случайная (III.A),2-Пропорциональная рулетка (III.B)\n- ";
        cin >> select_choice;
        cin.ignore( 10000,'\n');
        if (select_choice == 1 || select_choice==2) break;
        cout<<"неверный ввод."<<endl;
   } while (true);

    int mut_choice=1;
    do {
        cout<<"оператор мутации: 1-Простая (V.A),2-Инверсия (V.F),3-Случайный выбор из обоих\n- ";
        cin >> mut_choice;
        cin.ignore(10000,'\n');
        if (mut_choice>=1 && mut_choice <= 3) break;
        cout<<"неверный ввод."<<endl;
   } while (true);

    vector<int> rezlt;

    for (int gen=0; gen<Const::kolPok; gen++)
    {
        novPop.clear();

        vector<pair<int,int>> rodit;
        if (select_choice == 1)rodit = selSlu();
        else rodit=selShk( tekPop);

        if (cross_choice == 1) kroOdn(rodit , tekPop , novPop);
        else if (cross_choice == 2) kroDvu( rodit,tekPop ,novPop);
        else kroCik(rodit,tekPop,novPop);

        applyMutation(novPop,mut_choice);

        Osob lchPot=selElt(novPop);

        if (lchPot.znach != 0) {
            rezlt.push_back(lchPot.znach);

            int indPlh = 0;

            for (int i = 1;i < (int)tekPop.size(); i++) {
                if (tekPop[i].znach < tekPop[indPlh].znach)indPlh = i;
           }

            if (lchPot.znach > tekPop[indPlh].znach)tekPop[indPlh]=lchPot;
       }

   }

    int indLch=0;
    for (int i=1; i<(int)tekPop.size(); i++)
        if (tekPop[i].znach > tekPop[indLch].znach)
            indLch = i;
    if (!rezlt.empty())
        rezlt[ rezlt.size()-1 ]=tekPop[indLch].znach;

    int best_x_val = decode(tekPop[indLch].hrom);
    int best_f_val= tekPop[indLch].znach;
    cout<<"\nлучшее найденное знач: x="<<best_x_val<<",f(x)="<<best_f_val<<endl;
    cout<<"аналитический максимум f(10)="<<Func(10)<<endl;

    grfGr(rezlt);
    return 0;
}