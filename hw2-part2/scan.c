#include "scan.h"

bin_operator g_plus = NULL;
bin_operator g_cross = NULL;
bin_operator g_companion = NULL;
int gb_call_op_point = 0;

typedef struct Pair {
	double first;
	double second;
} Pair;

Pair op_point(const Pair z1, const Pair z2) {
	Pair ret;
	
	ret.first = g_companion(z1.first, z2.first);
	ret.second = g_plus(g_cross(z1.second, z2.first), z2.second);
	
	return ret;
}

int up_sweep(Pair* c, int size) {
	int d_max = (int) (log(size)/log(2));
	
	_Cilk_for (int d = 0; d < d_max; ++d) {
		int p = 1 << d; // pow(2,d)
		int i = 0;
		
	  _Cilk_for (int i = p-1; i < (size-p); i += p*2) {
			if (gb_call_op_point) {
				c[i+p] = op_point(c[i],c[i+p]);
			} else {
				c[i+p].first = g_plus(c[i].first,c[i+p].first);
			}
		}
	}
	
	return EXIT_SUCCESS;
	/* def __up_sweep(self):
		for d in range(int(math.log(len(self.c),2))):
			p = int(math.pow(2,d))
			i = p-1
			while i+p < len(self.c):
				self.c[i+p] = self.__op_point(self.c[i],self.c[i+p])
				i += p*2 */
}

int down_sweep(Pair* c, int size) {
	int max_d = log(size)/log(2) - 1;
	
  _Cilk_for (int d = max_d; d > 0; --d) {
		int p = 1 << (d-1);
		for (int i = (1 << d)-1; i < (size-p); i += p*2) {
			if (gb_call_op_point) {
				c[i+p] = op_point(c[i],c[i+p]);
			} else {
				c[i+p].first = g_plus(c[i].first,c[i+p].first);
			}
		}
	}
	
	return EXIT_SUCCESS;
	/* 
	def __down_sweep(self):
		for d in reversed(range(1,int(math.log(len(self.c),2)))):
			p = int(math.pow(2,d))
			i = p-1
			p /= 2
			while i+p < len(self.c):
				self.c[i+p] = self.__op_point(self.c[i],self.c[i+p])
				i += p*2
	 */
}

int scan(double* out, const double* a, const double* b, const unsigned int size, bin_operator plus, bin_operator cross, bin_operator companion) {
	Pair* c = NULL;
	int t_real_size = size;
  int t_size = 1 << (int)ceil(log(size)/log(2)); // make sure our array has a power of 2 len

	if ((c = (Pair*) malloc (sizeof(Pair)*t_size)) == NULL) {
		return EXIT_FAILURE;
	}
	
	gb_call_op_point = (cross != NULL);

  _Cilk_for (int i = 0; i < t_size; ++i) {
    c[i].first = i < t_real_size ? a[i] : 0;
    c[i].second = (gb_call_op_point && (i < t_real_size))? b[i] : 0;
  }
	
	g_plus = plus;
	g_cross = cross;
	g_companion = (companion == NULL) ? cross : companion;
	
	up_sweep(c, t_size);
	down_sweep(c, t_size);
	
	_Cilk_for (int i = 0; i < t_real_size; ++i) {
		out[i] = gb_call_op_point ? c[i].second : c[i].first;
	}
	
	free(c);
	return EXIT_SUCCESS;
}


