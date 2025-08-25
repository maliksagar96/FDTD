#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " file1 file2\n";
        return 1;
    }

    ifstream f1(argv[1]), f2(argv[2]);
    if (!f1.is_open() || !f2.is_open()) {
        cerr << "Error opening files.\n";
        return 1;
    }

    string line1, line2;
    int line_no = 0;
    bool same = true;

    while (true) {
        bool r1 = static_cast<bool>(getline(f1, line1));
        bool r2 = static_cast<bool>(getline(f2, line2));
        line_no++;

        if (!r1 && !r2) break; // both ended
        if (r1 != r2 || line1 != line2) {
            cout << "Difference at line " << line_no << endl;
            same = false;
            break; // stop at first difference
        }
    }

    if (same) cout << "Files are identical." << endl;
    else cout << "Files differ." << endl;

    return 0;
}
