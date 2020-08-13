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
        

if __name__=="__main__":
    unittest.main()
