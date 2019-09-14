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
# include "compress.h"

# define _CRT_SECURE_NO_WARNINGS

using namespace std;
const int MAX_CHA_NUM = 1 << 28;//maximum length of a chromosome
const int min_rep_len = 15;   //minimum replace length, matched string length exceeds min_rep_len, saved as matched information
const int LINE_CHA_NUM = 200; //maximum character number of one line
const int VEC_SIZE = 1 << 20; //length for other character arrays
const char integerDecoding[4] = { 'A', 'C', 'G', 'T' };//decoding 0,1,2,3 to A,C,G,T respectively
int *lineWidthArray;          //the array of line width, stores all the line widths of to-be-compressed sequences
int line_width_array_len = 0; //the length of the line width array

inline void initial()   // allocate memory
{
	ref_code = new char[MAX_CHA_NUM];
	seq_code = new char[MAX_CHA_NUM];
	ref_low = new POSITION_RANGE[VEC_SIZE];
	seq_low = new POSITION_RANGE[VEC_SIZE];
	diff_low = new POSITION_RANGE[VEC_SIZE];
	nCha = new POSITION_RANGE[VEC_SIZE];
	identifier_vec.reserve(seqNumber);
	matchResult.reserve(VEC_SIZE);
	matchResult_vec.reserve(seqNumber);
}

void decompressClear()   //release memory
{
	delete[] ref_code;
	delete[] seq_code;
	delete[] ref_low;
	delete[] seq_low;
	delete[] diff_low;
	delete[] low_loc;
	delete[] nCha;
	delete[] spe_cha;
	for (vector<char *>::iterator it = seqName.begin(); it != seqName.end(); it++)
		if (NULL != *it)
		{
			delete[] *it;
			*it = NULL;
		}
	seqName.clear();
}

//restore the lowercase character information
void seqLowercaseReading(int _seq_low_len, int _diff_low_len)
{
	int loc;
	for (int i = 0, j =0; i < _seq_low_len; i++)
	{
		loc = low_loc[i];
		if (loc == 0)
		{
			seq_low[i].begin = diff_low[j].begin;
			seq_low[i].length = diff_low[j++].length;
		}
		else
		{
			seq_low[i].begin = ref_low[loc].begin;
			seq_low[i].length = ref_low[loc].length;
		}
	}
}

void runLengthDecoding(FILE *fp, int *&vec, int &length, int tolerance)
{
	int code_len, temp;
	vector<int> code;
	fscanf(fp, "%d", &code_len);
	for (int i = 0; i < code_len; i++)
	{
		fscanf(fp, "%d", &temp);
		code.push_back(temp);
	}

	length = 0;
	for (int i = 1; i < code_len; i += 2)
		length += code[i];
	if (length > 0)
	{
		vec = new int[length];
		int k = 0;
		for (int i = 0; i < code_len; i += 2)
			for (int j = 0; j < code[i + 1]; j++)
				vec[k++] = code[i]+j*tolerance;
	}
}

void readIdentifierData(FILE *fp, vector<string> &vec)
{
	int code_len, temp_int;
	vector<int> code;
	char *temp_str = new char[LINE_CHA_NUM];
	fscanf(fp, "%d", &code_len);
	for (int i = 0; i < code_len; i++)
	{
		fscanf(fp, "%d", &temp_int);
		code.push_back(temp_int);
	}
	for (int i = 0; i < code_len; i++)
	{
		while (fgets(temp_str, LINE_CHA_NUM, fp) != NULL)
		{
			if (temp_str[0] == '\n') continue;
			identifier = temp_str;
			for (int j = 0; j < code[i]; j++)
				vec.push_back(identifier);
			break;
		}
	}
	delete temp_str;
}

void readPositionRangeData(FILE *fp, int &_vec_len, POSITION_RANGE *&_vec)
{
	fscanf(fp, "%d", &_vec_len);
	for (int i = 0; i < _vec_len; i++)
		fscanf(fp, "%d%d", &_vec[i].begin, &_vec[i].length);
}

void readSpeChaData(FILE *fp, int &_vec_len, POSITION_SPE_CHA *&_vec)
{
		_vec = new POSITION_SPE_CHA[_vec_len];
		for (int i = 0; i < _vec_len; i++)
			fscanf(fp, "%d", &_vec[i].pos);

		//read special character information
		vector<int> arr;
		int size, temp;
		fscanf(fp, "%d", &size);
		for (int i = 0; i < size; i++)
		{
			fscanf(fp, "%d", &temp);
			arr.push_back(temp);
		}
		if (size != 1)
		{
			int bit_num = ceil(log(size) / log(2));//the bit number of representing a special character
			int v_num = floor(32.0 / bit_num);//the number of characters can be represented in 4 bytes

			for (int i = 0; i < _vec_len; )
			{
				unsigned int v;
				fscanf(fp, "%u", &v);
				vector<int> temp_arr;

				int temp_i = i;
				for (int j = 0; j < v_num && temp_i < _vec_len; j++, temp_i++)
				{
					int mod = v % (1 << bit_num);
					v >>= bit_num;
					temp_arr.push_back(arr[mod]);
				}
				for (int j = temp_arr.size() - 1; j >= 0 && i < _vec_len; j--, i++)
					_vec[i].ch = temp_arr[j];
			}
		}
		else
		{
			for (int i = 0; i < _vec_len; i++)
				_vec[i].ch = arr[0];
		}
}

void readOtherData(FILE *fp, int &_seq_low_len, int &_nCha_len, int &_spe_cha_len)
{
	//read lowercase character information
	int flag;
	fscanf(fp, "%d", &flag);
	if (!flag)
		readPositionRangeData(fp, _seq_low_len, seq_low);
	else
	{
		runLengthDecoding(fp, low_loc, _seq_low_len, 1);
		readPositionRangeData(fp, diff_low_len, diff_low);
		seqLowercaseReading(_seq_low_len, diff_low_len);
	}
	//read n character information
	readPositionRangeData(fp, _nCha_len, nCha);

	//read special character information
	fscanf(fp, "%d", &_spe_cha_len);
	if (_spe_cha_len > 0)
	    readSpeChaData(fp, _spe_cha_len, spe_cha);
}

//get the space number of one line
int spaceNumber(char *str)
{
	int num = 0;
	for (int i = 0; str[i]; i++)
		if (str[i] == ' ')
			num++;
	return num;
}
void readFirstMatchResult(FILE *fp, vector <MatchEntry> &_mr)
{
	int _pos, _length;
	char temp_str[10240];
	string _misStr = "";
	MatchEntry _me;
	_mr.clear();

	fgets(temp_str, 10240, fp);//Throw out the line break between other data and matchentity
	fgets(temp_str, 10240, fp);
	while (temp_str[0] != '\n')
	{
		if (spaceNumber(temp_str))//there is space in this line, means it is matched result
		{
			sscanf(temp_str, "%d%d", &_pos, &_length);
			_me.pos = _pos;
			_me.length = _length;
			_me.misStr = _misStr;
			_mr.push_back(_me);
			_misStr = "";
		}
		else   //there is no space in this line, means it is mismatched string
			_misStr += temp_str;
		fgets(temp_str, 10240, fp);
	}
}
void getMatchResult(int _seq_id, int _pos, int _length, vector <MatchEntry> &_mr)
{
	for (int i = 0; i < _length; i++)
		_mr.push_back(matchResult_vec[_seq_id][_pos++]);
}
void readSecondMatchResult(FILE *fp, vector <MatchEntry> &_mr)
{
	int pre_seq_id = 0, _seq_id, pre_pos = 0, _pos, _length, spaceNum;
	char temp_str[10240];
	string _misStr = "";
	MatchEntry _me;
	_mr.clear();

	fgets(temp_str, 10240, fp);//Throw out the line break between other data and matchentity
	fgets(temp_str, 10240, fp);
	while (temp_str[0] != '\n')
	{
		spaceNum = spaceNumber(temp_str);
		if (spaceNum == 2)//the space number is 2, means it is the matched result of the second match
		{
			sscanf(temp_str, "%d%d%d", &_seq_id, &_pos, &_length);
			_seq_id += pre_seq_id;
			pre_seq_id = _seq_id;
			_pos += pre_pos;
			_length += 2;
			pre_pos = _pos + _length;
			getMatchResult(_seq_id, _pos, _length, _mr);
		}
		else if (spaceNum == 1)//the space number is 1, means it is the matched result of the first match
		{
			sscanf(temp_str, "%d%d", &_pos, &_length);
			_me.pos = _pos;
			_me.length = _length;
			_me.misStr = _misStr;
			_mr.push_back(_me);
			_misStr = "";
		}
		else   //the space number is 0, means it is the mismatched string of the first match
			_misStr += temp_str;
		fgets(temp_str, 10240, fp);
	}
}

//restore to-be-compressed sequence base code according to the matched result of the first match and the second match
void readTargetSequenceCode(vector <MatchEntry> &_mr, char *seq_code, int &seq_code_len)
{
	int _pos, pre_pos=0, cur_pos,  _length, _seq_code_len=0, str_len;
	const char *_misStr;
	for (unsigned int i =0; i < _mr.size();i++)
	{
		_pos = _mr[i].pos;
		cur_pos = _pos + pre_pos;
		_length = _mr[i].length + min_rep_len;
		pre_pos = cur_pos + _length;
		_misStr = _mr[i].misStr.c_str();
		//record misStr first ,and then record matchentity when compressing, so restore misStr first, and then matchentity when decompressing
		str_len = strlen(_misStr);
		for (int k = 0; k < str_len - 1; k++)
			seq_code[_seq_code_len++] = integerDecoding[_misStr[k] - '0'];
		for (int m = cur_pos, n = 0; n < _length; n++, m++)
			seq_code[_seq_code_len++] = ref_code[m];
	}
	seq_code[_seq_code_len] = '\0';
	seq_code_len = _seq_code_len;
}

//restore the to-be-compressed sequence file according to the base code, special character information, N character information and lowercase information
void saveSequenceFile(FILE *fp, int _seqnum, int &seq_code_len, int &_seq_low_len, int &_nCha_len, int &_spe_cha_len)
 {
	for (int i = 1; i < _spe_cha_len; i++)
		spe_cha[i].pos += spe_cha[i-1].pos;

	char *temp_seq = new char[MAX_CHA_NUM];

	strcpy(temp_seq, seq_code);
	int tt = 0, j = 0;

	for (int i = 0; i < _spe_cha_len; i++)
	{
		while (tt < spe_cha[i].pos && tt < seq_code_len)
		{
			seq_code[j++] = temp_seq[tt++];
		}
		seq_code[j++] = spe_cha[i].ch + 'A';
	}
	while (tt < seq_code_len) {
		seq_code[j++] = temp_seq[tt++];
	}
	seq_code[j] = '\0';
	seq_code_len = j;

	int str_len = 0;
	int r = 0;
	char *str = new char[MAX_CHA_NUM];

	for (int i = 0; i < _nCha_len; i++)
	{
		for (int j = 0; j < nCha[i].begin; j++)
			str[str_len++] = seq_code[r++];
		for (int j = 0; j < nCha[i].length; j++)
			str[str_len++] = 'N';
	}
	while (r < seq_code_len)
		str[str_len++] = seq_code[r++];
	str[str_len] = '\0';
	const char *_identifier = identifier_vec[_seqnum].c_str();
	fprintf(fp, "%s", _identifier);
	int k = 0;
	for (int i = 0; i < _seq_low_len; i++)
	{
		k += seq_low[i].begin;
		int temp = seq_low[i].length;
		for (int j = 0; j < temp; j++)
		{
			str[k] = tolower(str[k]);
			k++;
		}
	}
	int _lineWidth = lineWidthArray[_seqnum];
	int temp_seq_len = 0;
	for (int i = 0, m = 0; i < str_len; i++, m++)
	{
		if (m == _lineWidth)
		{
			temp_seq[temp_seq_len++] = '\n';
			m=0;
		}
		temp_seq[temp_seq_len++] = str[i];
	}
	temp_seq[temp_seq_len++] = '\n';
	temp_seq[temp_seq_len] = '\0';
	fprintf(fp, "%s", temp_seq);
	delete[] temp_seq;
	delete[] str;
}

void decompress(char *filename)
{
	printf("Info: Decompressing...Please wait for a moment.\n");
	sec_seq_num = ceil(percent * seqNumber / 100);
	initial();
	char temp_filename[100], resultFilename1[100], resultFilename2[100], targetFilename[100], cmd[100];
	extractFileName(filename, temp_filename);
	sprintf(resultFilename1, "%s.hrcm", temp_filename);
	sprintf(resultFilename2, "%s.desc", temp_filename);
	sprintf(cmd, "./7za x %s.7z", temp_filename);
	system(cmd);

	FILE *fp2 = fopen(resultFilename2, "r");
	if (NULL == fp2)
	{
		printf("Error: fail to open %s.\n", resultFilename2);
		exit(-1);
	}
	runLengthDecoding(fp2, lineWidthArray, line_width_array_len, 0);
	readIdentifierData(fp2, identifier_vec);
	fclose(fp2);

	referenceSequenceExtraction(seqName[0]);//reference sequence information extraction

	FILE *fp1 = fopen(resultFilename1, "r");
	if (NULL == fp1)
	{
		printf("Error: fail to open %s.\n", resultFilename1);
		exit(-1);
	}
	for (int i = 1; i < seqNumber; i++)
	{
		readOtherData(fp1, seq_low_len, nCha_len, spe_cha_len);
		if (i == 1) readFirstMatchResult(fp1, matchResult);
		else        readSecondMatchResult(fp1, matchResult);
		if (i <= sec_seq_num && i != seqNumber - 1) matchResult_vec.push_back(matchResult);
		readTargetSequenceCode(matchResult, seq_code, seq_code_len);

		extractFileName(seqName[i], temp_filename);
		sprintf(targetFilename, "%s.fasta", temp_filename);
		FILE *fp = fopen(targetFilename, "w");
		if (NULL == fp)
		{
			printf("Error: fail to open target sequence file %s. No. %d\n", targetFilename, i);
			exit(-1);
		}
		saveSequenceFile(fp, i - 1, seq_code_len, seq_low_len, nCha_len, spe_cha_len);
		fclose(fp);
		printf("Decompressed sequence %s (%d/%d).\n", targetFilename, i, seqNumber - 1);
	}
	fclose(fp1);
	sprintf(cmd, "rm -f %s %s", resultFilename1, resultFilename2);
	system(cmd); //for linux
	decompressClear();
}