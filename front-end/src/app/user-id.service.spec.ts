import { HttpClientTestingModule } from '@angular/common/http/testing';
import { TestBed } from '@angular/core/testing';

import { UserIdService } from './user-id.service';
import { first, firstValueFrom } from 'rxjs';
import { ClientService } from './client.service';

describe('UserIdService', () => {
  let service: UserIdService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [ HttpClientTestingModule ]
    });
    service = TestBed.inject(UserIdService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('setUserId', () => {
    it('should call client with local storage userId if it exists', async () => {
      spyOn(localStorage, 'getItem').and.returnValue('mock_user_id');
      const clientService: ClientService = TestBed.inject(ClientService);
      const clientServiceSpy = spyOn(clientService, 'getUserId')
      await service.setUserId();
      expect(clientServiceSpy).toHaveBeenCalledOnceWith('mock_user_id');
    });

    it('should call client with empty string if local storage userId does not exist', async () => {
      spyOn(localStorage, 'getItem').and.returnValue(null);
      const clientService: ClientService = TestBed.inject(ClientService);
      const clientServiceSpy = spyOn(clientService, 'getUserId')
      await service.setUserId();
      expect(clientServiceSpy).toHaveBeenCalledOnceWith('');
    });
    
    it('should update userId in local storage and in service as expected', async () => {
      const clientService: ClientService = TestBed.inject(ClientService);
      spyOn(clientService, 'getUserId').and.returnValue(Promise.resolve('mock_user_id'));
      const localStorageSpy = spyOn(localStorage, 'setItem');
      await service.setUserId();
      expect(localStorageSpy).toHaveBeenCalledOnceWith('userId', 'mock_user_id');
      const response = await firstValueFrom(service.userId.pipe(first()));
      expect(response).toEqual('mock_user_id');
    });
  });
});
