import { HttpTestingController, provideHttpClientTesting } from '@angular/common/http/testing';
import { MatSnackBar, MatSnackBarModule } from '@angular/material/snack-bar';
import { TestBed, fakeAsync, tick } from '@angular/core/testing';
import { provideHttpClient, withInterceptorsFromDi } from '@angular/common/http';

import { BackendHttpClient } from './backend-http.client';

describe('BackendHttpClient', () => {
  let service: BackendHttpClient;
  let snackBar: MatSnackBar;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        MatSnackBarModule
      ],
      providers: [
        provideHttpClient(withInterceptorsFromDi()),
        provideHttpClientTesting()
      ]
    });
    service = TestBed.inject(BackendHttpClient);
    snackBar = TestBed.inject(MatSnackBar);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('getUserId', () => {
    it('should request to the correct address and return the response', fakeAsync(async () => {
      const mockHttp = TestBed.inject(HttpTestingController);
      const return_value = service.getUserId('input_id');
      tick();
      const req = mockHttp.expectOne('http://localhost:5000/api/getuserid?user_id=input_id');
      expect(req.request.method).toBe('GET');
      req.flush({user_id: 'output_id'});
      mockHttp.verify();
      expect(await return_value).toBe('output_id');
    }));

    it('should open snack bar if the response is an error', fakeAsync(async () => {
      const mockHttp = TestBed.inject(HttpTestingController);
      const spy = spyOn(snackBar, 'open');
      service.getUserId('').then(async () => {
        expect(spy).toHaveBeenCalledOnceWith('Error retrieving User ID. Please refresh the page and try again.', 'Close', { duration: 5000 });
      });
      tick();
      const req = mockHttp.expectOne('http://localhost:5000/api/getuserid?user_id=');
      expect(req.request.method).toBe('GET');
      req.flush('', { status: 406, statusText: 'Bad Request'});
      mockHttp.verify();
    }));
  });

  describe('getDatasetInfo', () => {
    it('should request to the correct address and return the response', fakeAsync(async () => {
      const mockHttp = TestBed.inject(HttpTestingController);
      const return_value = service.getDatasetInfo();
      tick();
      const req = mockHttp.expectOne('http://localhost:5000/api/getdatasetinfo');
      expect(req.request.method).toBe('GET');
      req.flush({datasets: [
        {key: 'datasetKey', title: 'datasetText', samples: [], integration_obs: [], grouping_obs: []}
      ]});
      mockHttp.verify();
      expect(await return_value).toEqual([{key: 'datasetKey', title: 'datasetText', samples: [], integration_obs: [], grouping_obs: []}]);
    }));

    it('should open snack bar if the response is an error', fakeAsync(async () => {
      const mockHttp = TestBed.inject(HttpTestingController);
      const spy = spyOn(snackBar, 'open');
      service.getDatasetInfo().then(async () => {
        expect(spy).toHaveBeenCalledOnceWith('Error retrieving datasets. Please refresh the page and try again.', 'Close', { duration: 5000 });
      });
      tick();
      const req = mockHttp.expectOne('http://localhost:5000/api/getdatasetinfo');
      expect(req.request.method).toBe('GET');
      req.flush('', { status: 406, statusText: 'Bad Request'});
      mockHttp.verify();
    }));
  });
});
