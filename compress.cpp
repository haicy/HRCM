# include <iostream>
//# include <sys/timeb.h>  //windows
# include <sys/time.h>  //Linux
# include <cstring>
# include <vector>
//# include <windows.h>    //windows
# include <unistd.h>   //linux
# include <stdio.h>
# include <stdlib.h>
# include <cmath>
# include "decompress.h"
# define _CRT_SECURE_NO_WARNINGS

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
}MatchEntry;                  //for first match result



const int MAX_SEQ_NUM = 2000;//maximum sequence number
const int MAX_CHA_NUM = 1 << 28;//maximum length of a chromosome
const int LINE_CHA_NUM =200;    //maximum character number of one line
const int kMerLen = 14; //the length of k-mer
const int kmer_bit_num = 2 * kMerLen; //bit numbers of k-mer
const int hashTableLen = 1 << kmer_bit_num; // length of hash table
const int VEC_SIZE = 1 <<20; //length for other character arrays
const int min_rep_len = 15;   //minimum replace length, matched string length exceeds min_rep_len, saved as matched information

string identifier;
int lineWidth, ref_code_len, seq_code_len, ref_low_len, seq_low_len, diff_low_len, nCha_len, spe_cha_len, seqNumber, seqBucketLen;
int percent, sec_seq_num; //the percentage and reference sequence number used for the second compress 
char *ref_code, *seq_code;
char *dismatched_str; //mismatched subsequence
int *refLoc; //reference hash location
int *refBucket; //reference hash bucket
int *seqBucket;        //sequence hash bucket
int *low_loc;  //lowercase tuple location

POSITION_RANGE *ref_low, *seq_low, *diff_low, *nCha;
POSITION_SPE_CHA *spe_cha;


vector<char *> seqName;
vector <string> identifier_vec;
vector <int> lineWidth_vec;
vector <POSITION_RANGE> ref_nCha, seq_nCha, lowCha; //N character vector and lowercase character vector;
vector <MatchEntry> matchResult;                   //store the first match result of one sequence
vector <MatchEntry> misMatchEntry;                //store the mismatched match entity of the second match
vector < vector <MatchEntry> > matchResult_vec;   //store the match result of all sequences of the first match
vector <int *> seqBucket_vec;                    //store the hash bucket of all sequences
vector < vector<int> > seqLoc_vec;               //store the collision elements of all sequences when creating hash table


inline void initial()   // allocate memory
{
	ref_code = new char[MAX_CHA_NUM];
	seq_code = new char[MAX_CHA_NUM];
	refBucket = new int[hashTableLen];
	refLoc = new int[MAX_CHA_NUM];
	ref_low = new POSITION_RANGE[VEC_SIZE];
	seq_low = new POSITION_RANGE[VEC_SIZE];
	diff_low = new POSITION_RANGE[VEC_SIZE];
	low_loc = new int[VEC_SIZE];
	nCha = new POSITION_RANGE[VEC_SIZE];
	spe_cha = new POSITION_SPE_CHA[VEC_SIZE];
	identifier_vec.reserve(seqNumber);
	lineWidth_vec.reserve(seqNumber);
	seqBucket_vec.reserve(seqNumber);
	seqLoc_vec.reserve(seqNumber);
}

int integerCoding(char ch) { //encoding ACGT
	switch (ch) {
	case 'A': return 0;
	case 'C': return 1;
	case 'G': return 2;
	case 'T': return 3;
	default : return 4;
	}
}

//read the input file, get the filename of every sequence and the sequence amount
int readFile(char *filename)
{
	FILE* fp = fopen(filename, "r");
	if (NULL == fp) {
		printf("Error: failed to open file %s\n", filename);
		exit(-1);
	}
	char *temp_name = new char[LINE_CHA_NUM];  //the length of filename
	while (fscanf(fp,"%s", temp_name) != EOF)
	{
		seqName.push_back(temp_name);
		temp_name = new char[LINE_CHA_NUM];
	}
	seqNumber = seqName.size();
	delete temp_name;
	return seqNumber;
}

void referenceSequenceExtraction(char *str_referenceName)
{
	int _seq_code_len = 0, _ref_low_len = 1, letters_len = 0;//record lowercase from 1, diff_lowercase_loc[i]=0 means mismatching
	char temp_cha;
	bool flag = true;
	char cha[LINE_CHA_NUM];      //the content of one line

	printf("The reference file is %s.\n", str_referenceName);
	FILE* fp = fopen(str_referenceName, "r");
	if (NULL == fp)
	{
		printf("Error: fail to open reference file %s.\n", str_referenceName);
		exit(-1);
	}

	fgets(cha, LINE_CHA_NUM, fp);
	while (fscanf(fp, "%s", cha) != EOF)
	{
		for (unsigned int i = 0; i < strlen(cha); i++) {
			temp_cha = cha[i];
			if (islower(temp_cha))
			{
				if (flag) //previous is upper case
				{
					flag = false; //change status of flag
					ref_low[_ref_low_len].begin = letters_len;
					letters_len = 0;
				}
				temp_cha = toupper(temp_cha);
			}
			else            //this case is upper case
			{
				if (!flag)  //previous is lower case
				{
					flag = true;
					ref_low[_ref_low_len++].length = letters_len;
					letters_len = 0;
				}
			}
			if (temp_cha == 'A' || temp_cha == 'C' || temp_cha == 'G' || temp_cha == 'T')
				ref_code[_seq_code_len++] = temp_cha;
			letters_len++;
		}
	}
	if (!flag)  //if flag=false, don't forget record the length
		ref_low[_ref_low_len++].length = letters_len;

	fclose(fp);
	ref_code_len = _seq_code_len;
    ref_low_len = _ref_low_len - 1;
}

void targetSequenceExtraction(char *str_sequenceName, int &_seq_code_len, int &_seq_low_len, int &_nCha_len, int &_spe_cha_len)
{
	FILE *fp = fopen(str_sequenceName, "r");
	if (NULL == fp) {
		printf("Error: fail to open sequence file %s.\n", str_sequenceName);
		return;
	}

	_seq_code_len = _seq_low_len = _nCha_len = _spe_cha_len = 0;
	int letters_len = 0, n_letters_len = 0;
	bool flag = true, n_flag = false;
	char cha[LINE_CHA_NUM];      //the content of one line
	char temp_cha;

	//get the identifier
	fgets(cha, LINE_CHA_NUM, fp);
	identifier = cha;
	identifier_vec.push_back(identifier);

	//get the lineWidth
	if (fscanf(fp, "%s", cha) != EOF)
		lineWidth = strlen(cha);
	lineWidth_vec.push_back(lineWidth);
	fseek(fp, -1L * (lineWidth+1), SEEK_CUR);//set the 'fp' to the beginning of the file again

	while (fscanf(fp, "%s", cha) != EOF)
	{
		for (unsigned int i = 0; i < strlen(cha); i++)
		{
			temp_cha = cha[i];
			if (islower(temp_cha))
			{
				if (flag) //previous is upper case
				{
					flag = false;
					seq_low[_seq_low_len].begin = letters_len;
					letters_len = 0;
				}
				temp_cha = toupper(temp_cha);
			}
			else
			{
				if (isupper(temp_cha))
				{
					if (!flag)
					{
						flag = true;
						seq_low[_seq_low_len++].length = letters_len;
						letters_len = 0;
					}
				}
			}
			letters_len++;

			//temp_cha is an upper letter
			if (temp_cha == 'A' || temp_cha == 'C' || temp_cha == 'G' || temp_cha == 'T')
				seq_code[_seq_code_len++] = temp_cha;
			else if (temp_cha != 'N')
			{
				spe_cha[_spe_cha_len].pos = _seq_code_len;
				spe_cha[_spe_cha_len++].ch = temp_cha - 'A';
			}
			if (!n_flag)
			{
				if (temp_cha == 'N')
				{
					nCha[_nCha_len].begin = n_letters_len;
					n_letters_len = 0;
					n_flag = true;
				}
			}
			else
			{
				if (temp_cha != 'N')
				{
					nCha[_nCha_len++].length = n_letters_len;
					n_letters_len = 0;
					n_flag = false;
				}
			}
			n_letters_len++;
		}
	}

	if (!flag)
		seq_low[_seq_low_len++].length = letters_len;

	if (n_flag)
		nCha[_nCha_len++].length = n_letters_len;

	for (int i = _spe_cha_len - 1; i > 0; i--)
		spe_cha[i].pos -= spe_cha[i-1].pos;

	fclose(fp);
}

// construction of k-mer hashing table and linked list
void kMerHashingConstruct()//complete
{
	//initialize the point array
	for (int i = 0; i < hashTableLen; i++)
		refBucket[i] = -1;
	unsigned int value = 0;
	int step_len = ref_code_len - kMerLen + 1;

	//calculate the value of the first k-mer
	for (int k = kMerLen - 1; k >= 0; k--) {
		value <<= 2;
		value += integerCoding(ref_code[k]);
	}
	refLoc[0] = refBucket[value];
	refBucket[value] = 0;

	int shift_bit_num = (kMerLen * 2 - 2);
	int one_sub_str = kMerLen - 1;

	//calculate the value of the following k-mer using the last k-mer
	for (int i = 1; i < step_len; i++) {
		value >>= 2;
		value += (integerCoding(ref_code[i + one_sub_str])) << shift_bit_num;
		refLoc[i] = refBucket[value];    //refLoc[i] record the list of same values
		refBucket[value] = i;
	}
}

void codeFirstMatch(char *tar_seq_code, int tar_seq_len, vector<MatchEntry> &matchResult)
{
	//struct  timeval  c1_start;
	//struct  timeval  c1_end;
	//unsigned long c1_timer;
	//gettimeofday(&c1_start, NULL);
	int pre_pos = 0;
	int step_len = tar_seq_len - kMerLen + 1;
	int max_length, max_k;
	int i, id, k, ref_idx, tar_idx, length, cur_pos, tar_value;
	string mismatched_str;
	mismatched_str.reserve(10240);

	MatchEntry me;
	matchResult.reserve(VEC_SIZE);

	for (i = 0; i < step_len; i++) 
	{
		tar_value = 0;

		//calculate the hash value of the first k-mer
		for (k = kMerLen - 1; k >= 0; k--) 
		{
			tar_value <<= 2;
			tar_value += integerCoding(tar_seq_code[i + k]);
		}

		id = refBucket[tar_value];
		if (id > -1) 
		{                      //there is a same k-mer in ref_seq_code
			max_length = -1;
			max_k = -1;

			//search the longest match in the linked list
			for (k = id; k != -1; k = refLoc[k]) 
			{
					ref_idx = k + kMerLen;
					tar_idx = i + kMerLen;
					length = kMerLen;

					while (ref_idx < ref_code_len && tar_idx < tar_seq_len && ref_code[ref_idx++] == tar_seq_code[tar_idx++]) 
						length++;

					if (length >= min_rep_len && length > max_length)
					{
						max_length = length;
						max_k = k;
					}
			}

			if (max_length > -1) //exist a k-mer, its length is larger then min_rep_len
			{              
				//then save matched information
				cur_pos = max_k - pre_pos;      //delta-coding for cur_pos
				me.pos = cur_pos;
				me.length = max_length - min_rep_len;
				me.misStr = mismatched_str;
				matchResult.push_back(me);

				i += max_length;
				pre_pos = max_k + max_length;
				mismatched_str = "";
				if (i < tar_seq_len) 
					mismatched_str += '0' + integerCoding(tar_seq_code[i]);//mismatched_str stores the integer code of nucleotides
				continue;
			}
		}
		mismatched_str += '0' + integerCoding(tar_seq_code[i]);
	}
	if (i < tar_seq_len)
	{
		for (; i < tar_seq_len; i++)
			mismatched_str += '0' + integerCoding(tar_seq_code[i]);
		me.pos = 0;
		me.length = -min_rep_len;                //no match information, not 0 ,is -min_rep_len;
		me.misStr = mismatched_str;
		matchResult.push_back(me);
	}
	//gettimeofday(&c1_end, NULL);
	//c1_timer = 1000000 * (c1_end.tv_sec - c1_start.tv_sec) + c1_end.tv_usec - c1_start.tv_usec;
	//printf("codeFirstMatch time = %lf ms; %lf s\n", c1_timer / 1000.0, c1_timer / 1000.0 / 1000.0);
}

//lowercase character information matching
void seqLowercaseMatching(int _seq_low_len, int &_diff_low_len)
{
	int start_position = 1;
	_diff_low_len = 0;

	//initialize the diff_low_loc, diff_low_loc record the location of the same lowercase element
	memset(low_loc, 0, sizeof(int)*_seq_low_len);
	for (int i = 0; i < _seq_low_len; i++)
	{
	//search from the start_position to the end
	for (int j = start_position; j < ref_low_len; j++)
	{
		if ((seq_low[i].begin == ref_low[j].begin) && (seq_low[i].length == ref_low[j].length))
		{
			low_loc[i] = j;
			start_position = j + 1;
			break;
		}
	}

	//search from the start_position to the begin
	if (low_loc[i] == 0)
	{
		for (int j = start_position - 1; j > 0; j--) {
			if ((seq_low[i].begin == ref_low[j].begin) && (seq_low[i].length == ref_low[j].length))
			{
				low_loc[i] = j;
				start_position = j + 1;
				break;
			}
		}
	}

	//record the mismatched information
	if (low_loc[i] == 0)
	{
		diff_low[_diff_low_len].begin = seq_low[i].begin;
		diff_low[_diff_low_len++].length = seq_low[i].length;
	}
	}
}

void runLengthCoding(FILE *fp, int *vec, int length, int tolerance)
{
	vector<int> code;
	if (length > 0)
	{
		code.push_back(vec[0]);
		int cnt = 1;
		for (int i = 1; i < length; i++)
		{
			if (vec[i] - vec[i-1] == tolerance)
				cnt++;
			else
			{
				code.push_back(cnt);
				code.push_back(vec[i]);
				cnt = 1;
			}
		}
		code.push_back(cnt);
	}
	int code_len = code.size();
	fprintf(fp, "%d ", code_len);

	for (int i = 0; i < code_len; i++)
		fprintf(fp, "%d ", code[i]);
}

void saveIdentifierData(FILE *fp, vector<string> &vec)
{
	vector<string> _vec;
	vector<int> code;

	_vec.push_back(vec[0]);
	string pre_str = vec[0];
	int cnt = 1;
	for (unsigned int i = 1; i < vec.size(); i++)
	{
		if ((vec[i]) == pre_str)
			cnt++;
		else
		{
			code.push_back(cnt);
			_vec.push_back(vec[i]);
			pre_str = vec[i];
			cnt = 1;
		}
	}
	code.push_back(cnt);
	int code_len = code.size();
	fprintf(fp, " %d", code_len);
	for (int i = 0; i < code_len; i++)
		fprintf(fp, " %d", code[i]);
	fprintf(fp, "\n");
	for (int i = 0; i < code_len; i++)
		fprintf(fp, "%s", _vec[i].c_str());
}

void savePositionRangeData(FILE *fp, int _vec_len, POSITION_RANGE *_vec)
{
	fprintf(fp, "%d ", _vec_len);
	for (int i = 0; i < _vec_len; i++)
		fprintf(fp, "%d %d ", _vec[i].begin, _vec[i].length);
}

void saveSpeChaData(FILE *fp, int _vec_len, POSITION_SPE_CHA *_vec)
{
	int flag[26],temp;
	for (int i = 0; i < 26; i++)
		flag[i] = -1;
	vector<int> arr;
	for (int i = 0; i < _vec_len; i++)
	{
		fprintf(fp, "%d ", _vec[i].pos);
		temp = _vec[i].ch;
		if (flag[temp] == -1)
		{
			arr.push_back(temp);
			flag[temp] = arr.size() - 1;   //arr vector stores all types of special characters and flag array constructs the index
		}
	}

	int size = arr.size();       
	fprintf(fp, "%d ", size);             //size is the type amount of special characters
	for (int i = 0; i < size; i++)
		fprintf(fp, "%d ", arr[i]);
	if (size != 1)
	{
		unsigned int bit_num = ceil(log(size) / log(2));//the bit number of representing a special character
		unsigned int v_num = floor(32.0 / bit_num);     //the number of characters can be represented in 4 bytes
		for (int i = 0; i < _vec_len; )
		{
			unsigned int v = 0;
			for (unsigned int j = 0; j < v_num && i < _vec_len; j++, i++)
			{
				v <<= bit_num;
				v += flag[_vec[i].ch];
			}
			fprintf(fp, "%u ", v);
		}
	}
}

void saveOtherData(FILE *fp, int _seq_low_len, int _nCha_len, int _spe_cha_len)
{
	//save lowercase character information
	int flag = 0;
	if (_seq_low_len > 0 && ref_low_len > 0)
    {
        seqLowercaseMatching(_seq_low_len, diff_low_len);
        if (2 * diff_low_len < _seq_low_len)
        {
            flag = 1;
            fprintf(fp, "%d ", flag);
            runLengthCoding(fp, low_loc, _seq_low_len, 1);
            //printf("the length of diff_low_len is : %d ", diff_low_len);
            savePositionRangeData(fp, diff_low_len, diff_low);
        }
    }
    if(!flag)
    {
        fprintf(fp, "%d ", flag);
        savePositionRangeData(fp, _seq_low_len, seq_low);
    }

	//save n character information
    savePositionRangeData(fp, _nCha_len, nCha);

	//save special character information
	fprintf(fp, "%d ", _spe_cha_len);
	if (_spe_cha_len > 0)
	{
		saveSpeChaData(fp, _spe_cha_len, spe_cha);
	}
    fprintf(fp,"\n");//The end of other data
}

void saveMatchEntry(FILE *fp, MatchEntry &_me)
{
	if (!_me.misStr.empty())
		fprintf(fp, "%s\n", _me.misStr.c_str());
	fprintf(fp, "%d %d\n", _me.pos, _me.length);
}

void saveFirstMatchResult(FILE *fp, vector <MatchEntry> &_mr)
{
	for (unsigned int i = 0; i < _mr.size(); i++)
		saveMatchEntry(fp, _mr[i]);
    fprintf(fp, "\n");//The end flag of the first target sequence.
}

int getHashValue(MatchEntry &_me)
{
	int result = 0;
	for (unsigned int i = 0; i < _me.misStr.size(); i++)
		result += _me.misStr[i] * 92083;
	result += _me.pos * 69061 + _me.length * 51787;
	return result % seqBucketLen;
}

//get the nearest prime larger than number as the length of hash bucket
int getNextPrime(const int number)
{
	int cur = number + 1;
	bool prime = false;
	while (!prime)
	{
		prime = true;
		for (int l = 2; l < sqrt(number) + 1; l++)
		{
			if (cur%l == 0)
			{
				prime = false;
				break;
			}
		}

		if (!prime) cur++;
	}
	return cur;
}

//Construct the hash index for the second match
void matchResultHashConstruct(vector<MatchEntry> &_mr)
{
	int hashValue1, hashValue2, hashValue;
	vector<int> seqLoc;
	seqLoc.reserve(VEC_SIZE);
	int *seqBucket = new int[seqBucketLen];
	for (int i = 0; i < seqBucketLen; i++)
		seqBucket[i] = -1;

	//construct hash table
	hashValue1 = getHashValue(_mr[0]);
	if (_mr.size() < 2) hashValue2 = 0;
	else hashValue2 = getHashValue(_mr[1]);
	hashValue = abs(hashValue1 + hashValue2) % seqBucketLen;
	seqLoc.push_back(seqBucket[hashValue]);
	seqBucket[hashValue] = 0;
	
	for (unsigned int i = 1; i < _mr.size()-1; i++)
	{
		hashValue1 = hashValue2;
		hashValue2 = getHashValue(_mr[i+1]);
		hashValue = abs(hashValue1 + hashValue2) % seqBucketLen;
		seqLoc.push_back(seqBucket[hashValue]);
		seqBucket[hashValue] = i;
	}
	seqLoc_vec.push_back(seqLoc);
	seqBucket_vec.push_back(seqBucket);
	seqLoc.clear();
}

bool compareMatchEntry(MatchEntry &ref, MatchEntry &tar)//complete
{
	if (ref.pos == tar.pos && ref.length == tar.length && ref.misStr == tar.misStr)
		return true;
	else
		return false;
}

int getMatchLength(vector <MatchEntry> &ref_me, unsigned int ref_idx, vector <MatchEntry> &tar_me, unsigned int tar_idx)
{
	int length = 0;
	while (ref_idx < ref_me.size() && tar_idx < tar_me.size() && compareMatchEntry(ref_me[ref_idx++], tar_me[tar_idx++]))
		length++;
	return length;

}


//the second match
void codeSecondMatch(FILE *fp, vector<MatchEntry> &_mr, int seqNum)
{
	int hashValue;
	int pre_seq_id=1;
	int max_pos=0, pre_pos=0, delta_pos, length, max_length, delta_length, seq_id=0, delta_seq_id;
	int id, pos, secondMatchTotalLength=0;
	unsigned int i;

	for (i = 0; i < _mr.size()-1; i++)
	{
		if(_mr.size()<2) hashValue = abs(getHashValue(_mr[i])) % seqBucketLen;
		else hashValue = abs(getHashValue(_mr[i]) + getHashValue(_mr[i+1])) % seqBucketLen;
		max_length = 0;
		for (int m = 0; m < min( seqNum-1, sec_seq_num ); m++)
		{
			id = seqBucket_vec[m][hashValue];
			if (id!=-1)
			{
				for (pos = id; pos!=-1; pos = seqLoc_vec[m][pos])
				{
					length = getMatchLength(matchResult_vec[m], pos, _mr, i);
					if (length > 1 && length > max_length)
					{
						seq_id = m + 1;  
						max_pos = pos;
						max_length = length;
					}
				}
			}
		}
		if (max_length)
		{
			delta_seq_id = seq_id - pre_seq_id;        //delta encoding
			delta_length = max_length - 2;            //delta encoding, the minimum replace length of the second match is 2
			delta_pos = max_pos - pre_pos;            //delta encoding
			pre_seq_id = seq_id;
			pre_pos = max_pos + max_length;
			secondMatchTotalLength += max_length;

			//firstly save mismatched matchentry
			if (!misMatchEntry.empty())
			{
				for (unsigned int i = 0; i < misMatchEntry.size(); i++)
					saveMatchEntry(fp, misMatchEntry[i]);
				misMatchEntry.clear();   
			}

			//secondly save matched matchentry
			fprintf(fp, "%d %d %d\n", delta_seq_id, delta_pos, delta_length);
			i += max_length - 1;
		}
		else
		{
			misMatchEntry.push_back(_mr[i]);
		}
	}
	if (i == _mr.size()-1)  misMatchEntry.push_back(_mr[i]);
	if (!misMatchEntry.empty())
	{
		for (unsigned int j = 0; j < misMatchEntry.size(); j++)
			saveMatchEntry(fp, misMatchEntry[j]);
		misMatchEntry.clear();
	}
    fprintf(fp, "\n");//The end flag of target sequence.
}

void codeMatch(FILE *fp)
{
	vector<MatchEntry> mr;
	sec_seq_num = ceil(percent * seqNumber / 100);    //get the amount of reference sequences for the second match
	seqBucketLen = getNextPrime(VEC_SIZE);
	for (int i = 1; i < seqNumber; i++)
	{
		targetSequenceExtraction(seqName[i], seq_code_len, seq_low_len, nCha_len, spe_cha_len);
		saveOtherData(fp, seq_low_len, nCha_len, spe_cha_len);
		codeFirstMatch(seq_code, seq_code_len, mr);
		if (i <= sec_seq_num)
		{
			matchResult_vec.push_back(mr);
			matchResultHashConstruct(mr);
		}
		if (i == 1)                                   //save the match result of the first to-be-compressed sequence directly
			saveFirstMatchResult(fp, mr);
		else
			codeSecondMatch(fp, mr, i);
		printf("Compressed sequence %s (%d/%d).\n", seqName[i], i, seqNumber - 1);
		mr.clear();
	}
}
void compressClear()  	//release memory
{
	delete[] ref_code;
	delete[] seq_code;
	delete[] refLoc;
	delete[] refBucket;
	delete[] ref_low;
	delete[] seq_low;
	delete[] diff_low;
	delete[] low_loc;
	delete[] nCha;
	delete[] spe_cha;

	for (vector<int *>::iterator it = seqBucket_vec.begin(); it != seqBucket_vec.end(); it++)
		if (NULL != *it)
		{
			delete[] *it;
			*it = NULL;
		}
	seqBucket_vec.clear();
	for (vector<char *>::iterator it = seqName.begin(); it != seqName.end(); it++)
		if (NULL != *it)
		{
			delete[] *it;
			*it = NULL;
		}
	seqName.clear();
}

void extractFileName(char *srcFileName, char *destFileName)
{
	unsigned int i=0, j=0;
	for (; i < strlen(srcFileName); i++) //get the characters between '/' and '.'
	{
		if (srcFileName[i] == '/')
		{
			j = 0;
			continue;
		}
		else
		{
			if (srcFileName[i] == '.')
				break;
			else
				destFileName[j++] = srcFileName[i];
		}
	}
	destFileName[j] = '\0';
}

void compress(char *filename)
{
	
	if (seqNumber > 1)
	{
		printf("Info: Compressing...Please wait for a moment.\n");
		initial();
		char temp_filename[100], resultFilename1[100], resultFilename2[100], cmd[100];
		extractFileName(filename, temp_filename);
		sprintf(resultFilename1, "%s.hrcm", temp_filename);//stores the match result of all to-be-compressed sequences
		sprintf(resultFilename2, "%s.desc", temp_filename);//stores the identifier and line width of all to-be-compressed sequences

		referenceSequenceExtraction(seqName[0]);//reference sequence information extraction
		kMerHashingConstruct();//construct the hash index for the first match based on k-mer

		FILE *fp1 = fopen(resultFilename1, "w");
		if (NULL == fp1) {
			printf("Error: fail to open %s.\n", resultFilename1);
			exit(-1);
		}
		codeMatch(fp1);
		fclose(fp1);

		FILE *fp2 = fopen(resultFilename2, "w");
		if (NULL == fp2) {
			printf("Error: fail to open %s.\n", resultFilename2);
			exit(-1);
		}
		runLengthCoding(fp2, &lineWidth_vec[0], seqNumber - 1, 0);//save lineWidth data
		saveIdentifierData(fp2, identifier_vec);//save identifier data
		fclose(fp2);
		sprintf(cmd, "./7za a %s.7z %s %s -m0=PPMd", temp_filename, resultFilename1, resultFilename2);
		system(cmd); //for linux
		sprintf(cmd, "rm -f %s %s", resultFilename1, resultFilename2);
		system(cmd); //for linux
		compressClear();
	}
	else
		printf("Error: There is not any to-be-compressed sequence, nothing to be done.\n");
	
}

void show_usage() {
	cout << "HRCM v1.0\n";
	cout << "Usage: hrcm {compress | decompress}  -r {ref-file-path}{ [-t] {tar-file-path}|[-f] {filename} [percent]}\n";
	cout << "  {compress | decompress} is mode,  choose one of them according to requirement, required\n";
	cout << "  -r is the reference, the {ref-file-path} followed, required\n";
	cout << "  -t is the target, a single to-be-compressed file path {tar-file-path} followed, optional\n";
	cout << "  -f is the alternative option of -t, a set of to-be-compressed file paths included in {filename}, optional\n";
    cout << "  [percent] is the percentage of the second-level matching, default is 10, means 10 percent of sequences will be used for the second-level matching, optional when -f, illegal when -t\n";
	cout << "Examples:\n";
	cout << "  hrcm compress -r hg17_chr21.fa -t hg18_chr21.fa\n";
	cout << "  hrcm decompress -r hg17_chr21.fa -t hg18_chr21_ref_hg17_chr21.fa\n";
	cout << "  hrcm compress -r hg17_chr21.fa -f filename.txt 20\n";
	cout << "  hrcm decompress -r hg17_chr21.fa -f filename.txt 20\n";
}

int main(int argc, char *argv[]) {
	bool arg_flag = true, ref_flag = false, tar_flag = false, tar_set_flag = false, compressed = false, decompressed = false;
	char *mode = NULL, *ref_file = NULL, *tar_file = NULL, *per = NULL;
	int oc;

	struct  timeval  start;
	struct  timeval  end;
	unsigned long timer;
	gettimeofday(&start,NULL);//first argument get the time result, second argument get the timezone

	if (argc < 6 || argc > 7)
	{
		show_usage();
		return 0;
	}
	mode = argv[1];                         //mode = compress or decompress
	percent = 10; //default is 10
	if (argc == 7)
	{
		per = argv[6];
		percent = atoi(per); 
	}
	while ((oc = getopt(argc, argv, "r:t:f:")) >= 0)
	{
		switch (oc)
		{
		case 'r':
			ref_file = optarg;
			ref_flag = true;
			break;
		case 't':
			tar_file = optarg;
			tar_flag = true;
			break;
		case 'f':
			tar_file = optarg;
			tar_set_flag = true;
			break;
		case '?':
			arg_flag = false;
			break;
		}
	}
	if (!(tar_flag ^ tar_set_flag))
		arg_flag = false;
	if (arg_flag && ref_flag)
	{
		char *temp_name = new char[LINE_CHA_NUM];
		strcpy(temp_name, ref_file);
		seqName.push_back(temp_name);
		if (tar_flag)
		{
			temp_name = new char[LINE_CHA_NUM];
			strcpy(temp_name, tar_file);
			seqName.push_back(temp_name);
			seqNumber = 2;
		}
		else
			readFile(tar_file);		
		if (strcmp(mode, "compress") == 0)
		{
			compress(tar_file);
			compressed = true;
		}
		if (strcmp(mode, "decompress") == 0)
		{
			decompress(tar_file);
			decompressed = true;
		}

	}
	if (!(compressed||decompressed))
	{
		show_usage();
		return 0;
	}
	gettimeofday(&end,NULL);
	timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	if (compressed)
		printf("Total compression time = %lf ms; %lf s\n", timer/1000.0, timer/1000.0/1000.0);
	else
		printf("Total decompression time = %lf ms; %lf s\n", timer/1000.0, timer/1000.0/1000.0);
}

