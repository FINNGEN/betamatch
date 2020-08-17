pipeline {
  agent any

  stages {
    stage('Build') {
      steps {
        script {
          c = docker.build("phewas-development/betamatch:test-" + "$BUILD_NUMBER", "-f docker/Dockerfile ./")
          docker.withRegistry('http://gcr.io/phewas-development', 'gcr:phewas-development') {
            c.push("test-${env.BUILD_NUMBER}")
          }
    }
      }
    }
    stage('Tests') {
      /*set up tests*/
      steps{
        c = docker.build("phewas-development/betamatch:test-" + "$BUILD_NUMBER", "-f docker/Dockerfile ./")
        c.inside("-u root"){
          sh """
          python3 -m pip install pylint pytest safety pyflakes mypy prospector bandit
          #dl cromwell 48 womtool
          curl -O https://github.com/broadinstitute/cromwell/releases/download/48/womtool-48.jar
          chmod +x womtool-48.jar
          #run it
          ./womtool-48.jar betamatch/wdl/betamatch_github.wdl -i betamatch/wdl/betamatch_github.wdl
          #run python tests
          cd betamatch
          python3 -m pytest 
          """
        }
        /*sh 'python --version'
        sh 'python3 -m pip install pylint pytest safety pyflakes mypy prospector bandit'
        sh 'curl -O https://github.com/broadinstitute/cromwell/releases/download/48/womtool-48.jar'
        sh 'chmod +x womtool-48.jar'
        sh './womtool-48.jar betamatch/wdl/betamatch_github.wdl -i betamatch/wdl/betamatch_github.wdl'
        sh 'cd betamatch'
        sh 'python3 -m pytest'*/
      }
    }
  /*
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
  */
  }
}