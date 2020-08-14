pipeline {
  agent any

  stages {
    stage('Build') {
      steps {
	script {  c = docker.build("betamatch/jenkins_tests:test-" + "$BUILD_NUMBER", "-f docker/Dockerfile ./")
	  docker.withRegistry('http://eu.gcr.io/phewas-development', 'gcr:phewas-development') {
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