#ifndef WELCOME_H
#define WELCOME_H

#include <iostream>
#include <string>

#define VERSION "0.0.1"
#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "?"
#endif

using namespace std;

const string compile_date = __DATE__;
const string compile_time = __TIME__;
const string git_commit = GIT_COMMIT_HASH;


void printWelcomeCoulombIntegral(int rank=0)
{
if (rank!=0) return;
const string banner =  "\
\n\
  ______                _                _        _____                                           _ \n\
 / _____)              | |              | |      (_____)        _                                | |\n\
| /        ___   _   _ | |  ___   ____  | | _       _    ____  | |_    ____   ____   ____   ____ | |\n\
| |       / _ \\ | | | || | / _ \\ |    \\ | || \\     | |  |  _ \\ |  _)  / _  ) / _  | / ___) / _  || |\n\
| \\_____ | |_| || |_| || || |_| || | | || |_) )   _| |_ | | | || |__ ( (/ / ( ( | || |    ( ( | || |\n\
 \\______) \\___/  \\____||_| \\___/ |_|_|_||____/   (_____)|_| |_| \\___) \\____) \\_|| ||_|     \\_||_||_|\n\
                                                                            (_____|                 \n\
";

    cout << "\ncompiled from git commit:  " << git_commit << " on " << compile_date << " at " << compile_time << "\n";
    cout << banner;
    //cout << "                                                                                         version "<< VERSION << "\n\n\n";
    cout << "\nCitation:\tMerkel, K. & Ortmann, F. Journal of Physics: Materials 7, 015001 (2023)\n";
    cout << "         \thttps://dx.doi.org/10.1088/2515-7639/ad06cd\n\n\n";
}

/*void printWelcomeOpticalDipole(int rank)
{
if (rank!=0) return;
const string banner =  "\
\n\
_______          _____ _____                ______       ________ _____                 ______      \n\
__  __ \\________ __  /____(_)_____________ ____  /       ___  __ \\___(_)________ ______ ___  /_____ \n\
_  / / /___  __ \\_  __/__  / _  ___/_  __ `/__  /        __  / / /__  / ___  __ \\_  __ \\__  / _  _ \\\n\
/ /_/ / __  /_/ // /_  _  /  / /__  / /_/ / _  /         _  /_/ / _  /  __  /_/ // /_/ /_  /  /  __/\n\
\\____/  _  .___/ \\__/  /_/   \\___/  \\__,_/  /_/          /_____/  /_/   _  .___/ \\____/ /_/   \\___/ \n\
        /_/                                                             /_/                         \n\
";
    cout << "\ncompiled from git commit:  " << git_commit << " on " << compile_date << " at " << compile_time << "\n";
    cout << banner;
    cout << "                                                                                         version "<< VERSION << "\n\n\n";
}
*/
#endif // WELCOME_H