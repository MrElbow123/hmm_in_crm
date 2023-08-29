#include <Rcpp.h>
using namespace Rcpp;

// This is a c++ version of the the forward part of the forward backward algorithm for a categorical emission distribution
// Where forward is extracted from Baum-Welch algorithm


//' @keywords internal
// [[Rcpp::export(rng = false)]]
List cat_mult_fw_cpp(NumericMatrix allprobs, NumericMatrix gamma, int m, int n, NumericVector delta) {

  int i, t, k;
  NumericVector foo(m); // 创建一个m位的向量foo
  NumericVector foo1(m);
  double sumfoo;
  NumericMatrix alpha_prob(m,n); //创建一个m*n的向量alpha_prob——2*900
  double lscale;
  NumericMatrix lalpha(m, n);

  // 1 alpha初始化 在alpha_prob中----------

  foo = delta * allprobs(0,_); //delta * allprobs第0行 两个2*1的矩阵相乘
  sumfoo = 0;
  for (i = 0; i < m; i++){
    sumfoo += foo[i];
  }
  for (i = 0; i < m; i++){
    alpha_prob(i,0) = foo[i]/sumfoo;
  }
  lscale = log(sumfoo);
  for(i = 0; i < m; i++){
    lalpha(i, 0) = log(alpha_prob(i, 0)) + lscale;
  }


  // 2 alpha迭代----------

  for (t = 1; t < n; t++){
    for (i = 0; i < m; i++){
      foo1[i] = 0;
      for (k = 0; k < m; k++){
        foo1[i] += alpha_prob(k, (t - 1)) * gamma(k,i); // alpha*A
      }
    }
    foo = foo1 * allprobs(t,_); // B*alpha*A
    sumfoo = 0;
    for (i = 0; i < m; i++){
      sumfoo += foo[i];
    }
    for (i = 0; i < m; i++){
      alpha_prob(i,t) = foo[i]/sumfoo;
    }
    lscale = lscale + log(sumfoo);
    for(i = 0; i < m; i++){
      lalpha(i, t) = log(alpha_prob(i, t)) + lscale;
    }
  }

  return List::create(alpha_prob, lalpha);
  // alpha(i,t)表示给定参数组合delta/gamma/beta的情况下，观测序列是y1,y2..,yt，且t时刻状态为i的概率
  // alpha_prob是2*900矩阵 通过forward[[1]][,t] 取得t时刻的alpha(1,t)和alpha(2,t) 取最终的就是t=900
  // lalpha是做了对数处理
  }
