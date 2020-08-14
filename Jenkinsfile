pipeline {
  agent any

  stages {
    stage('Build') {
      steps {
	script {  c = docker.build("phewas-development/betamatch:test-" + "$BUILD_NUMBER", "-f docker/Dockerfile ./")
	  docker.withRegistry('http://gcr.io/phewas-development', 'gcr:phewas-development') {
			      c.push("test-${env.BUILD_NUMBER}")
			      }
	}
      }
    }
    stage('Metrics') {
      steps {
	withSonarQubeEnv('sonar') {

	  sh "${tool("sonar")}/bin/sonar-scanner \
          -Dsonar.projectKey=${JOB_NAME} \
          -Dsonar.sources=. \
          -Dsonar.css.node=. \
          -Dsonar.host.url=${DEFAULT_SONAR_URL} \
          -Dsonar.login=${SONAR_LOGIN}"

	}
      }
    }
  }
}