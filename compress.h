#ifndef COMPRESS_H_INCLUDED
#define COMPRESS_H_INCLUDED
# include <cstring>
# include <vector>
using namespace std;

typedef struct
{
	int begin;
	int length;
} POSITION_RANGE;             //for N character and lowercase character

typedef struct
{
	int pos;
	int ch;
} POSITION_SPE_CHA;         //for special character

typedef struct {
	int pos;
	int length;
	string misStr;
}MatchEntry;                  //for the first match result

extern string identifier;
extern int ref_code_len, seq_code_len, ref_low_len, seq_low_len, diff_low_len, nCha_len, spe_cha_len, seqNumber;
extern int percent, sec_seq_num; //the percentage and referenced sequence number used for second compress 
extern char *dismatched_str, *ref_code, *seq_code;
extern POSITION_RANGE *ref_low, *seq_low, *diff_low, *nCha;
extern POSITION_SPE_CHA *spe_cha;
extern int *low_loc; //lowercase tuple location

extern vector <char *> seqName;
extern vector <string>  identifier_vec;
extern vector <POSITION_RANGE> ref_nCha, seq_nCha, lowCha; //N character vector and lowercase character vector;
extern vector <MatchEntry> matchResult;                 //store the first match result of one sequence
extern vector <MatchEntry> misMatchEntry;               //store the mismatched match entity of the second match
extern vector < vector <MatchEntry> > matchResult_vec;  //store the match result of all sequences of the first match

void extractFileName(char *srcFileName, char *destFileName);
void referenceSequenceExtraction(char *str_referenceName);
int readFile(char *filename); 

#endif // COMPRESS_H_INCLUDED
