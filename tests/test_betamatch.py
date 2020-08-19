import unittest
import betamatch

class TestBetamatch(unittest.TestCase):
    def test_cov(self):
        """Test weighted covariance
        Cases:
            Unweighted covariance (population)
            Weighted covariance
                TODO
        """
        #unweighted covariance
        vec1 = [0.47668557, 0.41799971, 0.58757653, 0.86832309, 0.17516729, 0.92833054, 0.09948352, 0.05089571, 0.67195539, 0.04547867]
        vec2 = [0.06116168, 0.63266968, 0.86774159, 0.08985116, 0.50030305, 0.96073755, 0.84801626, 0.94017402, 0.21918019, 0.56216323]
        weights = [1.0]*10
        cov = betamatch.weighted_cov(vec1,vec2,weights)
        validate = -0.02896266 # calculated elsewhere. Note that R's cov calculates sample covariance, which gives slightly different answer.
        self.assertAlmostEqual(cov,validate)

    def test_r2(self):
        """Test correlation calculation
        """
        #y [1] -0.8111950  0.9137861  1.3570089  3.1348695  5.6801601  4.2758522  7.9149642  7.9932773  8.4275412  9.5862034 10.7976482
        #z [1]  0.3414767  0.3372171  0.6494706  4.5005896  6.2254523  4.9553089  7.4915341  8.5150125  7.9974688  8.6979310 12.0419944
        #cor (y,z) = 0.9760983
        validation=0.9760983 #from R using cor
        y = [-0.8111950, 0.9137861, 1.3570089, 3.1348695, 5.6801601, 4.2758522, 7.9149642, 7.9932773, 8.4275412, 9.5862034, 10.7976482]
        z = [0.3414767, 0.3372171, 0.6494706, 4.5005896, 6.2254523, 4.9553089, 7.4915341, 8.5150125, 7.9974688, 8.6979310, 12.0419944]
        w = [1.0]*11
        cor = betamatch.weighted_pearsonr(y,z,w)
        self.assertAlmostEqual(cor,validation,places=6)

if __name__=="__main__":
    unittest.main()
