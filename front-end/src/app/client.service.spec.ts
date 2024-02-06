import { HttpClientTestingModule, HttpTestingController } from '@angular/common/http/testing';
import { TestBed, fakeAsync, tick } from '@angular/core/testing';

import { ClientService } from './client.service';

describe('ClientService', () => {
  let service: ClientService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        HttpClientTestingModule,
      ]
    });
    service = TestBed.inject(ClientService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('getUserId', () => {
    it('should request to the correct address and return the response', fakeAsync(async () => {
      const mockHttp = TestBed.inject(HttpTestingController);
      const return_value = service.getUserId('input_id');
      tick();
      const req = mockHttp.expectOne('http://127.0.0.1:5000/getuserid?user_id=input_id');
      expect(req.request.method).toBe('GET');
      req.flush({text: 'output_id'});
      mockHttp.verify();
      expect(await return_value).toBe('output_id');
    }));
  });
});
